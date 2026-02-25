from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import DataLoader, Dataset, random_split
from tqdm import tqdm

from p53cad.core.logging import get_logger
from p53cad.core.runtime import select_device
from p53cad.engine.latent import ManifoldEmbedder


class DMSDataset(Dataset):
    def __init__(
        self,
        sequences: List[str],
        scores: List[float],
        embedder: ManifoldEmbedder,
        per_position: bool = False,
    ):
        self.sequences = sequences
        self.scores = torch.tensor(scores, dtype=torch.float32)
        self.embedder = embedder
        self.per_position = per_position
        self.cache: Dict[str, torch.Tensor] = {}

    def __len__(self):
        return len(self.sequences)

    def __getitem__(self, idx):
        seq = self.sequences[idx]
        score = self.scores[idx]

        if seq in self.cache:
            emb = self.cache[seq]
        else:
            with torch.no_grad():
                raw = self.embedder.encode(seq)  # (1, L, D)
                if self.per_position:
                    emb = raw.squeeze(0).cpu()  # (L, D)
                else:
                    emb = raw.mean(dim=1).squeeze(0).cpu()  # (D,)
            self.cache[seq] = emb

        return emb, score


class FunctionalNet(nn.Module):
    """
    Single-head MLP oracle.
    """

    def __init__(
        self,
        input_dim: int = 1280,
        hidden_dim: int = 128,
        num_layers: int = 2,
        dropout: float = 0.2,
    ):
        super().__init__()
        depth = max(int(num_layers), 1)
        width = max(int(hidden_dim), 16)
        drop = float(max(0.0, min(dropout, 0.9)))

        # Preserve legacy default geometry for checkpoint compatibility:
        # hidden_dim=128, num_layers=2 -> [128, 64]
        hidden_dims = [width]
        for layer_idx in range(1, depth):
            next_width = max(int(round(width / (2 ** layer_idx))), 16)
            hidden_dims.append(next_width)

        layers: list[nn.Module] = []
        in_dim = input_dim
        for idx, out_dim in enumerate(hidden_dims):
            layers.append(nn.Linear(in_dim, out_dim))
            layers.append(nn.ReLU())
            # Match legacy behavior: dropout only between hidden layers.
            if idx < len(hidden_dims) - 1 and drop > 0:
                layers.append(nn.Dropout(drop))
            in_dim = out_dim
        layers.append(nn.Linear(in_dim, 1))
        self.net = nn.Sequential(*layers)

    def forward(self, x):
        target_device = self.net[0].weight.device
        if x.device != target_device:
            x = x.to(target_device)
        return self.net(x)


class AttentionPoolingNet(nn.Module):
    """Oracle with learnable attention pooling over sequence positions.

    Accepts (batch, seq_len, D) and learns which positions matter for
    functional scoring, eliminating the mean-pooling bottleneck.
    Also accepts (batch, D) for backward compatibility (passes through MLP head).

    Uses delta encoding: subtracts a cached WT baseline from per-position
    inputs so that mutated positions stand out as nonzero residuals.
    Without this, single-AA mutations are invisible among 393 identical positions.
    """

    def __init__(
        self,
        input_dim: int = 1280,
        hidden_dim: int = 256,
        n_heads: int = 4,
        dropout: float = 0.2,
    ):
        super().__init__()
        self.input_dim = input_dim

        # WT baseline for delta encoding (set during training, saved in state_dict)
        self.register_buffer("wt_baseline", None)

        # Learnable query vector for attention pooling
        self.query = nn.Parameter(torch.randn(1, 1, input_dim) * 0.02)

        # Multi-head attention: query attends over sequence positions
        self.attn = nn.MultiheadAttention(
            embed_dim=input_dim, num_heads=n_heads, dropout=dropout, batch_first=True
        )
        self.layer_norm = nn.LayerNorm(input_dim)

        # MLP head: input_dim → hidden → hidden//2 → 1
        self.head = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.ReLU(),
            nn.Linear(hidden_dim // 2, 1),
        )

    def set_wt_baseline(self, wt_emb: torch.Tensor) -> None:
        """Register WT per-position embedding for delta encoding.

        Parameters
        ----------
        wt_emb : Tensor
            Shape (L, D) or (1, L, D) — WT hidden states from ESM-2.
        """
        if wt_emb.dim() == 3:
            wt_emb = wt_emb.squeeze(0)
        self.wt_baseline = wt_emb.detach().clone()

    def forward(self, x):
        target_device = self.query.device
        if x.device != target_device:
            x = x.to(target_device)

        # Backward compat: if already pooled (batch, D), just run MLP head
        if x.dim() == 2:
            return self.head(x)

        # Delta encoding: subtract WT baseline so mutation positions have nonzero signal
        if self.wt_baseline is not None:
            wt = self.wt_baseline.to(x.device)
            x = x - wt.unsqueeze(0)  # (batch, L, D) - (1, L, D)

        # x: (batch, seq_len, D)
        batch_size = x.size(0)
        query = self.query.expand(batch_size, -1, -1)  # (batch, 1, D)

        # Attention pooling: query attends to all positions
        pooled, _attn_weights = self.attn(query, x, x)  # (batch, 1, D)
        pooled = self.layer_norm(pooled.squeeze(1))  # (batch, D)

        return self.head(pooled)


class FunctionalOracle:
    def __init__(
        self,
        model_path: Path | str = None,
        input_dim: int = 1280,
        hidden_dim: int = 128,
        num_layers: int = 2,
        dropout: float = 0.2,
        use_rtl: bool = False,
        arch: str = "legacy_mlp",
        device: Optional[str] = None,
    ):
        self.logger = get_logger(__name__)
        self.device = select_device(device)
        self.input_dim = input_dim
        self.hidden_dim = hidden_dim
        self.num_layers = max(int(num_layers), 1)
        self.dropout = float(max(0.0, min(dropout, 0.9)))
        self.arch_name = arch

        if arch == "attention_pooling":
            self.model: nn.Module = AttentionPoolingNet(
                input_dim=self.input_dim,
                hidden_dim=self.hidden_dim,
                n_heads=4,
                dropout=self.dropout,
            )
        else:
            self.model: nn.Module = FunctionalNet(
                input_dim=self.input_dim,
                hidden_dim=self.hidden_dim,
                num_layers=self.num_layers,
                dropout=self.dropout,
            )

        if use_rtl:
            self.logger.warning("RTL oracle was removed. Falling back to %s.", self.arch_name)

        if model_path and Path(model_path).exists():
            self._load_model(model_path)

        self.model = self.model.to(self.device)
        self.model.eval()

        # torch.compile requires Triton (Linux-only); skip on Windows/MPS
        import platform as _platform
        if self.device.type == "cuda" and _platform.system() != "Windows":
            try:
                import torch as _torch
                self.model = _torch.compile(self.model, mode="reduce-overhead")
                self.logger.info("Oracle compiled with torch.compile (reduce-overhead)")
            except Exception as exc:
                self.logger.warning("torch.compile skipped for oracle: %s", exc)

    def _load_model(self, model_path: Path | str) -> None:
        self.logger.info(f"Loading oracle from {model_path}")
        payload = torch.load(model_path, map_location=self.device)

        # New checkpoint format with explicit metadata.
        if isinstance(payload, dict) and "model_state_dict" in payload:
            arch = payload.get("arch", "legacy_mlp")
            self.input_dim = int(payload.get("input_dim", self.input_dim))
            self.hidden_dim = int(payload.get("hidden_dim", self.hidden_dim))
            self.num_layers = int(payload.get("num_layers", self.num_layers))
            self.dropout = float(payload.get("dropout", self.dropout))

            if arch == "attention_pooling":
                model = AttentionPoolingNet(
                    input_dim=self.input_dim,
                    hidden_dim=self.hidden_dim,
                    n_heads=int(payload.get("n_heads", 4)),
                    dropout=self.dropout,
                )
            elif arch == "legacy_mlp":
                model = FunctionalNet(
                    input_dim=self.input_dim,
                    hidden_dim=self.hidden_dim,
                    num_layers=self.num_layers,
                    dropout=self.dropout,
                )
            else:
                self.logger.warning(
                    "Unsupported oracle architecture '%s' in %s. "
                    "Falling back to legacy_mlp.",
                    arch,
                    model_path,
                )
                model = FunctionalNet(
                    input_dim=self.input_dim,
                    hidden_dim=self.hidden_dim,
                    num_layers=self.num_layers,
                    dropout=self.dropout,
                )
                self.model = model
                self.arch_name = "legacy_mlp"
                return

            try:
                sd = payload["model_state_dict"]
                # wt_baseline registered as None isn't in the expected state_dict;
                # pre-register with matching shape so load_state_dict can fill it.
                if "wt_baseline" in sd and hasattr(model, "wt_baseline") and model.wt_baseline is None:
                    model.register_buffer("wt_baseline", torch.zeros_like(sd["wt_baseline"]))
                model.load_state_dict(sd)
            except Exception as exc:
                self.logger.warning(
                    "Failed to load oracle checkpoint %s: %s. Using fresh %s weights.",
                    model_path,
                    exc,
                    arch,
                )
            else:
                self.model = model
                self.arch_name = arch
                self.logger.info(
                    "Loaded oracle architecture: %s (input=%d, hidden=%d, layers=%d, dropout=%.2f)",
                    arch,
                    self.input_dim,
                    self.hidden_dim,
                    self.num_layers,
                    self.dropout,
                )
            return

        # Legacy checkpoint format: raw state_dict only.
        state_dict = payload
        model = FunctionalNet(
            input_dim=self.input_dim,
            hidden_dim=self.hidden_dim,
            num_layers=self.num_layers,
            dropout=self.dropout,
        )
        try:
            model.load_state_dict(state_dict)
        except Exception as exc:
            self.logger.warning(
                "Failed to load legacy oracle checkpoint %s: %s. Using fresh legacy_mlp weights.",
                model_path,
                exc,
            )
            return

        self.model = model
        self.arch_name = "legacy_mlp"
        self.logger.info("Loaded legacy oracle checkpoint.")

    def train(
        self,
        dms_data,
        embedder: ManifoldEmbedder,
        epochs: int = 10,
        save_path: Path = None,
        val_split: float = 0.1,
        early_stopping_patience: int = 8,
        min_delta: float = 1e-4,
        batch_size: int = 32,
        seed: int = 42,
    ):
        self.logger.info(f"Training FunctionalOracle on {len(dms_data)} sequences...")
        self.model.train()

        use_per_position = isinstance(self.model, AttentionPoolingNet)

        # For attention oracle: compute WT baseline for delta encoding
        if use_per_position:
            from p53cad.data.dms import P53_WT

            with torch.no_grad():
                wt_emb = embedder.encode(P53_WT).squeeze(0).cpu()  # (L, D)
            self.model.set_wt_baseline(wt_emb)
            self.logger.info(
                "Set WT baseline for delta encoding (%s)", tuple(wt_emb.shape)
            )

        dataset = DMSDataset(
            dms_data["sequence"].tolist(),
            dms_data["score"].tolist(),
            embedder,
            per_position=use_per_position,
        )
        total_samples = len(dataset)
        val_fraction = float(max(0.0, min(val_split, 0.49)))
        patience = max(int(early_stopping_patience), 0)
        delta = float(max(min_delta, 0.0))
        train_loader: DataLoader
        val_loader: Optional[DataLoader] = None

        if total_samples >= 5 and val_fraction > 0.0:
            val_size = int(round(total_samples * val_fraction))
            val_size = max(1, min(val_size, total_samples - 1))
            train_size = total_samples - val_size
            split_gen = torch.Generator().manual_seed(seed)
            train_ds, val_ds = random_split(dataset, [train_size, val_size], generator=split_gen)
            train_loader = DataLoader(train_ds, batch_size=batch_size, shuffle=True)
            val_loader = DataLoader(val_ds, batch_size=batch_size, shuffle=False)
            self.logger.info(
                "Using validation split: train=%d, val=%d (%.1f%%)",
                train_size,
                val_size,
                100.0 * val_size / total_samples,
            )
        else:
            train_loader = DataLoader(dataset, batch_size=batch_size, shuffle=True)
            self.logger.info(
                "Validation disabled (samples=%d, val_split=%.3f). Training on full set.",
                total_samples,
                val_fraction,
            )

        optimizer = optim.AdamW(self.model.parameters(), lr=1e-3)
        criterion = nn.MSELoss()
        best_metric = float("inf")
        best_epoch = 0
        epochs_completed = 0
        patience_counter = 0
        best_state_dict: Optional[Dict[str, torch.Tensor]] = None

        for epoch in range(epochs):
            total_loss = 0.0

            for X, y in tqdm(train_loader, desc=f"Epoch {epoch + 1}"):
                X = X.to(self.device)
                y = y.to(self.device).unsqueeze(1)

                optimizer.zero_grad()
                pred = self.model(X)
                loss = criterion(pred, y)
                loss.backward()
                optimizer.step()
                total_loss += float(loss.item())

            epochs_completed = epoch + 1
            n_train_batches = max(len(train_loader), 1)
            train_loss = total_loss / n_train_batches
            metric_for_selection = train_loss
            val_loss: Optional[float] = None

            if val_loader is not None:
                self.model.eval()
                val_total = 0.0
                with torch.no_grad():
                    for X_val, y_val in val_loader:
                        X_val = X_val.to(self.device)
                        y_val = y_val.to(self.device).unsqueeze(1)
                        pred_val = self.model(X_val)
                        val_total += float(criterion(pred_val, y_val).item())
                self.model.train()
                n_val_batches = max(len(val_loader), 1)
                val_loss = val_total / n_val_batches
                metric_for_selection = val_loss
                self.logger.info(
                    "Epoch %d Train Loss: %.4f | Val Loss: %.4f",
                    epoch + 1,
                    train_loss,
                    val_loss,
                )
            else:
                self.logger.info("Epoch %d Train Loss: %.4f", epoch + 1, train_loss)

            if metric_for_selection + delta < best_metric:
                best_metric = metric_for_selection
                best_epoch = epoch + 1
                best_state_dict = {
                    name: tensor.detach().cpu().clone()
                    for name, tensor in self.model.state_dict().items()
                }
                patience_counter = 0
            elif val_loader is not None and patience > 0:
                patience_counter += 1
                if patience_counter >= patience:
                    self.logger.info(
                        "Early stopping at epoch %d (best epoch %d, best val loss %.4f).",
                        epoch + 1,
                        best_epoch,
                        best_metric,
                    )
                    break

        if best_state_dict is not None:
            self.model.load_state_dict(best_state_dict)
            self.logger.info(
                "Restored best checkpoint from epoch %d (selection metric %.4f).",
                best_epoch,
                best_metric,
            )
        self.model.eval()
        if save_path:
            checkpoint = {
                "arch": self.arch_name,
                "input_dim": self.input_dim,
                "hidden_dim": int(self.hidden_dim),
                "num_layers": int(self.num_layers),
                "dropout": float(self.dropout),
                "n_heads": int(self.model.attn.num_heads) if hasattr(self.model, "attn") else 0,
                "model_state_dict": self.model.state_dict(),
                "training": {
                    "epochs_requested": int(epochs),
                    "epochs_completed": int(epochs_completed),
                    "best_epoch": int(best_epoch) if best_epoch else int(epochs_completed),
                    "best_selection_metric": float(best_metric),
                    "val_split": float(val_fraction),
                    "early_stopping_patience": int(patience),
                    "min_delta": float(delta),
                    "train_samples": int(total_samples - (len(val_loader.dataset) if val_loader is not None else 0)),
                    "val_samples": int(len(val_loader.dataset) if val_loader is not None else 0),
                },
            }
            torch.save(checkpoint, save_path)
            self.logger.info(f"Saved oracle checkpoint to {save_path} (arch={self.arch_name})")

    def train_multimut(
        self,
        dms_data,
        embedder: ManifoldEmbedder,
        n_synthetic: int = 50_000,
        k_range: tuple[int, int] = (2, 15),
        epochs: int = 10,
        save_path=None,
        val_split: float = 0.1,
        early_stopping_patience: int = 5,
        batch_size: int = 32,
        lr: float = 2e-4,
        seed: int = 42,
    ):
        """Fine-tune oracle on synthetic multi-mutation augmented data.

        Generates n_synthetic additive pseudo-labeled combinations from the DMS
        single-mutation data, mixes them with the original singles (1:1 ratio),
        and fine-tunes the oracle at a reduced learning rate.  This exposes the
        attention mechanism to inputs with multiple simultaneously nonzero
        positions in the delta embedding — the regime encountered during
        inference on 18–30-mutation rescue candidates.

        Pseudo-labels use the thermodynamic additivity model (standard null model
        in double-mutant cycle analysis) with soft saturation.  Literature-
        validated rescues (N239Y+R249S, H168R+R249S, etc.) are injected as
        high-confidence anchors (repeated 5× for upweighting).

        Parameters
        ----------
        dms_data : pd.DataFrame
            Original single-mutation DMS data (sequences + scores).
        embedder : ManifoldEmbedder
            ESM-2 wrapper (must be the same model used for original training).
        n_synthetic : int
            Number of synthetic k-mutation combinations to generate.
        k_range : tuple[int, int]
            Range of mutations per synthetic example (uniform sampling).
        epochs : int
            Fine-tuning epochs.
        save_path : Path or str, optional
            Where to save the updated checkpoint.
        val_split : float
            Fraction of augmented dataset reserved for validation.
        early_stopping_patience : int
            Stop if val loss doesn't improve for this many epochs.
        batch_size : int
            Batch size for fine-tuning.
        lr : float
            Learning rate (lower than original training to prevent catastrophic
            forgetting of single-mutation knowledge).
        seed : int
            Random seed for synthetic data generation and splits.
        """
        import pandas as pd

        if not isinstance(self.model, AttentionPoolingNet):
            raise RuntimeError(
                "train_multimut requires AttentionPoolingNet oracle "
                "(--arch attention_pooling). Current arch: %s" % self.arch_name
            )

        self.logger.info(
            "Generating %d synthetic multi-mutation training pairs (k=%d–%d)...",
            n_synthetic, k_range[0], k_range[1],
        )
        syn_seqs, syn_scores = generate_multimut_augmentation(
            dms_data, n_samples=n_synthetic, k_range=k_range,
            seed=seed, include_literature=True,
        )
        self.logger.info("Generated %d synthetic pairs (including literature anchors).", len(syn_seqs))

        # Mix: original singles (full weight) + synthetic (full weight).
        # Using separate dataset objects so the DataLoader shuffles across both.
        orig_seqs = dms_data["sequence"].tolist()
        orig_scores = dms_data["score"].tolist()
        all_seqs = orig_seqs + syn_seqs
        all_scores = orig_scores + syn_scores
        self.logger.info(
            "Augmented dataset: %d original + %d synthetic = %d total.",
            len(orig_seqs), len(syn_seqs), len(all_seqs),
        )

        # Ensure WT baseline is set
        if self.model.wt_baseline is None:
            from p53cad.data.dms import P53_WT
            with torch.no_grad():
                wt_emb = embedder.encode(P53_WT).squeeze(0).cpu()
            self.model.set_wt_baseline(wt_emb)
            self.logger.info("Set WT baseline for delta encoding.")

        aug_df = pd.DataFrame({"sequence": all_seqs, "score": all_scores})

        # Re-use standard train() with lower LR and the augmented dataset.
        # Temporarily override the optimizer LR by monkey-patching a closure.
        original_adamw = optim.AdamW

        class _LRAdam(original_adamw):
            def __init__(self_inner, params, **kwargs):
                kwargs["lr"] = lr
                super().__init__(params, **kwargs)

        optim.AdamW = _LRAdam  # temporary swap
        try:
            self.train(
                dms_data=aug_df,
                embedder=embedder,
                epochs=epochs,
                save_path=save_path,
                val_split=val_split,
                early_stopping_patience=early_stopping_patience,
                batch_size=batch_size,
                seed=seed,
            )
        finally:
            optim.AdamW = original_adamw  # restore

        self.logger.info("Multi-mutation fine-tuning complete.")

    def _prepare_input(self, embedding: torch.Tensor) -> torch.Tensor:
        if not torch.is_tensor(embedding):
            embedding = torch.tensor(embedding, dtype=torch.float32)
        vec = embedding.to(self.device).float()
        if vec.dim() == 3:
            # Attention oracle consumes (batch, L, D) directly
            if isinstance(self.model, AttentionPoolingNet):
                return vec
            return vec.mean(dim=1)
        if vec.dim() == 2:
            return vec
        if vec.dim() == 1:
            return vec.unsqueeze(0)
        raise ValueError(f"Unsupported embedding shape for oracle prediction: {tuple(vec.shape)}")

    def predict(self, embedding: torch.Tensor) -> float:
        """Predict function score from an embedding tensor."""
        with torch.no_grad():
            vec = self._prepare_input(embedding)
            return float(self.model(vec).squeeze(-1).mean().item())

    def predict_batch(self, embeddings: torch.Tensor) -> torch.Tensor:
        """Predict function scores for a batch of embeddings."""
        vec = self._prepare_input(embeddings)
        return self.model(vec)

    def predict_with_routing(self, embedding: torch.Tensor) -> Dict[str, float]:
        """
        Backward-compatible shim: routing diagnostics are removed.
        """
        return {"score": self.predict(embedding), "arch": self.arch_name}


# ===========================================================================
# Multi-mutation training augmentation
# ===========================================================================

# Known literature multi-mutation rescues with experimentally estimated outcomes.
# Format: (mutations_list, z_score_estimate, citation)
# Z-score sign: positive = functional (oracle convention after DMS negation).
_LITERATURE_MULTIMUT: list[tuple[list[str], float, str]] = [
    # N239Y rescues R249S: restores L3 loop contacts (Nikolova et al. 2000)
    (["R249S", "N239Y"], 0.60, "Nikolova2000"),
    # H168R rescues R249S: mimics lost Arg249 sidechain (Baroni et al. 2004)
    (["R249S", "H168R"], 0.45, "Baroni2004"),
    # N268D global stabilizer: rescues R175H partial (Brachmann et al. 1998)
    (["R175H", "N268D"], 0.30, "Brachmann1998"),
    # T284R rescues R175H: restores zinc coordination geometry (Otsuka et al. 2007)
    (["R175H", "T284R"], 0.50, "Otsuka2007"),
    # N239Y global suppressor: also rescues R248Q (Baroni et al. 2004)
    (["R248Q", "N239Y"], 0.35, "Baroni2004"),
    # N239Y + N268D double stabilizer on R249S (Brachmann1998 + Nikolova2000)
    (["R249S", "N239Y", "N268D"], 0.75, "Nikolova2000+Brachmann1998"),
    # WT p53 reference: zero mutations = baseline (constructed)
    ([], 0.04, "WT_reference"),
]


def _apply_mutations_to_wt(mutations: list[str], wt: str) -> str:
    """Apply a list of mutation strings (e.g. 'R175H') to WT sequence."""
    from p53cad.data.dms import parse_single_mutation
    seq = list(wt)
    for m in mutations:
        parsed = parse_single_mutation(m)
        if parsed is not None:
            _, pos, mut_aa = parsed
            if 1 <= pos <= len(seq):
                seq[pos - 1] = mut_aa
    return "".join(seq)


def _additive_zscore(z_scores: list[float]) -> float:
    """Additive pseudo-label with soft saturation.

    Thermodynamic independence model: DMS Z-scores are additive (ΔΔG is additive
    for non-contacting mutations). Soft saturation prevents unrealistic values far
    outside the DMS training range (|Z| > 3.5).

    Formula: z_sum / (1 + |z_sum| / 7.0) keeps output in roughly (-4, +4).
    """
    if not z_scores:
        return 0.04  # WT baseline
    z_sum = float(sum(z_scores))
    return z_sum / (1.0 + abs(z_sum) / 7.0)


def generate_multimut_augmentation(
    dms_data,
    n_samples: int = 50_000,
    k_range: tuple[int, int] = (2, 15),
    seed: int = 42,
    include_literature: bool = True,
) -> tuple[list[str], list[float]]:
    """Generate synthetic multi-mutation training pairs via additive pseudo-labeling.

    Randomly samples k mutations from the DMS dataset (k uniform in k_range),
    combines their sequences, and computes an additive Z-score pseudo-label.
    Also injects known literature multi-mutation rescues as high-confidence anchors.

    The purpose is to expose the oracle's attention mechanism to inputs where
    multiple positions are nonzero in the delta embedding — the distributional
    regime it encounters during inference on 18–30-mutation rescue candidates.

    Parameters
    ----------
    dms_data : pd.DataFrame
        DMS dataset with columns 'sequence' and 'score'.
    n_samples : int
        Number of synthetic combinations to generate.
    k_range : tuple[int, int]
        Minimum and maximum number of mutations per synthetic example.
    seed : int
        Random seed for reproducibility.
    include_literature : bool
        Whether to add experimentally-validated multi-mutation anchors.

    Returns
    -------
    sequences : list[str]
    scores : list[float]
    """
    import random
    from p53cad.data.dms import P53_WT

    rng = random.Random(seed)

    # Build lookup: sequence → score (use all available DMS single-mutation data)
    seq_list = dms_data["sequence"].tolist()
    score_list = dms_data["score"].tolist()

    # Identify which positions are mutated in each DMS entry (vs WT)
    mut_positions: list[int] = []
    mut_aas: list[str] = []
    for seq in seq_list:
        changed = [i for i, (a, b) in enumerate(zip(P53_WT, seq)) if a != b]
        mut_positions.append(changed[0] if len(changed) == 1 else -1)
        mut_aas.append(seq[changed[0]] if len(changed) == 1 else "")

    # Only use confirmed single-mutation entries
    valid_indices = [i for i, p in enumerate(mut_positions) if p >= 0]

    synthetic_seqs: list[str] = []
    synthetic_scores: list[float] = []
    seen: set[str] = set()  # prevent duplicates

    k_lo, k_hi = int(k_range[0]), int(k_range[1])

    for _ in range(n_samples):
        k = rng.randint(k_lo, min(k_hi, len(valid_indices)))
        # Sample k DISTINCT positions to avoid conflicting mutations at same site
        sampled = rng.sample(valid_indices, k)
        positions_used: set[int] = set()
        accepted: list[int] = []
        for idx in sampled:
            p = mut_positions[idx]
            if p not in positions_used:
                positions_used.add(p)
                accepted.append(idx)
        if len(accepted) < 2:
            continue

        # Build combined sequence by applying all k mutations to WT
        seq_chars = list(P53_WT)
        z_parts: list[float] = []
        for idx in accepted:
            p = mut_positions[idx]
            seq_chars[p] = mut_aas[idx]
            z_parts.append(score_list[idx])

        combined_seq = "".join(seq_chars)
        if combined_seq in seen:
            continue
        seen.add(combined_seq)

        z_pseudo = _additive_zscore(z_parts)
        synthetic_seqs.append(combined_seq)
        synthetic_scores.append(z_pseudo)

    # Inject literature-validated multi-mutation anchors (repeated 5× for upweighting)
    if include_literature:
        for mutations, z_lit, _cite in _LITERATURE_MULTIMUT:
            if not mutations:
                lit_seq = P53_WT
            else:
                lit_seq = _apply_mutations_to_wt(mutations, P53_WT)
            for _ in range(5):
                synthetic_seqs.append(lit_seq)
                synthetic_scores.append(z_lit)

    return synthetic_seqs, synthetic_scores


def compute_masked_marginal_pll(
    embedder: ManifoldEmbedder,
    sequence: str,
    mutated_positions: List[int],
) -> float:
    """Compute masked marginal pseudo-log-likelihood at mutated positions.

    For each position in *mutated_positions* (1-indexed), mask the token,
    run ESM-2 forward, and accumulate log P(true_aa | context).  This is
    the gold-standard zero-shot fitness predictor (Meier et al. 2021).

    Too expensive for the inner optimisation loop (one forward pass per
    position) but ideal for ranking ~30 shortlist candidates.

    Parameters
    ----------
    embedder : ManifoldEmbedder
        Loaded ESM-2 wrapper with tokenizer and model.
    sequence : str
        Full-length protein sequence (e.g. 393 AA for p53).
    mutated_positions : list of int
        1-indexed residue positions to evaluate.

    Returns
    -------
    float
        Sum of log P(aa_i | context) over mutated positions.
    """
    if not mutated_positions:
        return 0.0

    # Tokenize without special tokens — positions map 1:1 to token indices
    inputs = embedder.tokenizer(
        sequence, return_tensors="pt", add_special_tokens=False
    ).to(embedder.device)
    input_ids = inputs.input_ids  # (1, L)

    mask_token_id = embedder.tokenizer.mask_token_id
    total_log_prob = 0.0

    with torch.no_grad():
        for pos in mutated_positions:
            idx = pos - 1  # 0-indexed
            if idx < 0 or idx >= input_ids.shape[1]:
                continue

            masked_ids = input_ids.clone()
            masked_ids[0, idx] = mask_token_id

            outputs = embedder.model(input_ids=masked_ids, return_dict=True)
            logits = outputs.logits  # (1, L, vocab)
            log_probs = F.log_softmax(logits[0, idx], dim=-1)
            true_token = input_ids[0, idx].item()
            total_log_prob += float(log_probs[true_token].item())

    return total_log_prob


def compute_conditional_rescue_scores(
    embedder: ManifoldEmbedder,
    cancer_seq: str,
    positions: List[int],
) -> Dict[int, Dict[str, float]]:
    """Compute ESM-2 P(aa | cancer_context) at each position.

    Batches all masked positions into a single GPU forward pass for efficiency.
    Uses the embedder's model directly to avoid dtype conflicts.

    Parameters
    ----------
    embedder : ManifoldEmbedder
        Loaded ESM-2 wrapper.
    cancer_seq : str
        Full-length cancer-mutant protein sequence.
    positions : list of int
        1-indexed residue positions to evaluate.

    Returns
    -------
    dict
        {pos: {aa_letter: log_prob, ...}, ...}
    """
    if not positions:
        return {}

    # Build batch of masked sequences (all masks at once)
    inputs = embedder.tokenizer(
        cancer_seq, return_tensors="pt", add_special_tokens=False
    ).to(embedder.device)
    input_ids = inputs.input_ids  # (1, L)
    mask_token_id = embedder.tokenizer.mask_token_id

    # Create batch: one sequence per position with that position masked
    batch_size = len(positions)
    if batch_size == 0:
        return {}
    
    # Build batched input - clone original and mask each position
    batch_input_ids = input_ids.repeat(batch_size, 1)  # (batch, L)
    for batch_idx, pos in enumerate(positions):
        idx = pos - 1
        if 0 <= idx < batch_input_ids.shape[1]:
            batch_input_ids[batch_idx, idx] = mask_token_id

    # Build AA token id → letter mapping
    aa_letters = "ACDEFGHIKLMNPQRSTVWY"
    aa_token_ids = []
    for aa in aa_letters:
        tid = embedder.tokenizer.convert_tokens_to_ids(aa)
        aa_token_ids.append(tid)

    result: Dict[int, Dict[str, float]] = {}

    # Single forward pass with all masked sequences
    with torch.no_grad():
        try:
            # Use embedder's model directly - this is the same way encode() works
            outputs = embedder.model.esm(batch_input_ids, output_hidden_states=True, return_dict=True)
            logits = embedder.model.lm_head(outputs.last_hidden_state)  # (batch, L, vocab)
            log_probs = torch.log_softmax(logits, dim=-1)
            
            # Extract scores for each masked position
            for batch_idx, pos in enumerate(positions):
                idx = pos - 1
                if 0 <= idx < log_probs.shape[1]:
                    pos_log_probs = log_probs[batch_idx, idx, :]  # (vocab,)
                    pos_scores: Dict[str, float] = {}
                    for aa_char, tid in zip(aa_letters, aa_token_ids):
                        if tid < pos_log_probs.shape[0]:
                            pos_scores[aa_char] = float(pos_log_probs[tid].item())
                    result[pos] = pos_scores
                    
        except Exception as e:
            # If GPU fails, return empty dict (campaign continues without this feature)
            import logging
            logging.getLogger(__name__).warning(f"Conditional rescue GPU compute failed: {e}")
            return {}

    return result
