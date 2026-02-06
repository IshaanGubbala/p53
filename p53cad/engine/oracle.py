from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional

import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, Dataset, random_split
from tqdm import tqdm

from p53cad.core.logging import get_logger
from p53cad.engine.latent import ManifoldEmbedder


class DMSDataset(Dataset):
    def __init__(self, sequences: List[str], scores: List[float], embedder: ManifoldEmbedder):
        self.sequences = sequences
        self.scores = torch.tensor(scores, dtype=torch.float32)
        self.embedder = embedder
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
                emb = raw.mean(dim=1).squeeze(0).cpu()  # (D,)
            self.cache[seq] = emb

        return emb, score


class FunctionalNet(nn.Module):
    """
    Single-head MLP oracle.
    """

    def __init__(
        self,
        input_dim: int = 320,
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


class FunctionalOracle:
    def __init__(
        self,
        model_path: Path | str = None,
        input_dim: int = 320,
        hidden_dim: int = 128,
        num_layers: int = 2,
        dropout: float = 0.2,
        use_rtl: bool = False,
    ):
        self.logger = get_logger(__name__)
        self.device = torch.device("mps" if torch.backends.mps.is_available() else "cpu")
        self.input_dim = input_dim
        self.hidden_dim = hidden_dim
        self.num_layers = max(int(num_layers), 1)
        self.dropout = float(max(0.0, min(dropout, 0.9)))
        self.arch_name = "legacy_mlp"
        self.model: nn.Module = FunctionalNet(
            input_dim=self.input_dim,
            hidden_dim=self.hidden_dim,
            num_layers=self.num_layers,
            dropout=self.dropout,
        )

        if use_rtl:
            self.logger.warning("RTL oracle was removed. Falling back to legacy_mlp.")

        if model_path and Path(model_path).exists():
            self._load_model(model_path)

        self.model = self.model.to(self.device)
        self.model.eval()

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
            model = FunctionalNet(
                input_dim=self.input_dim,
                hidden_dim=self.hidden_dim,
                num_layers=self.num_layers,
                dropout=self.dropout,
            )

            if arch != "legacy_mlp":
                self.logger.warning(
                    "Unsupported oracle architecture '%s' in %s. "
                    "RTL checkpoints are no longer used; retrain to produce legacy_mlp.",
                    arch,
                    model_path,
                )
                self.model = model
                self.arch_name = "legacy_mlp"
                return

            try:
                model.load_state_dict(payload["model_state_dict"])
            except Exception as exc:
                self.logger.warning(
                    "Failed to load oracle checkpoint %s: %s. Using fresh legacy_mlp weights.",
                    model_path,
                    exc,
                )
            else:
                self.model = model
                self.arch_name = "legacy_mlp"
                self.logger.info(
                    "Loaded oracle architecture: legacy_mlp (input=%d, hidden=%d, layers=%d, dropout=%.2f)",
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

        dataset = DMSDataset(
            dms_data["sequence"].tolist(),
            dms_data["score"].tolist(),
            embedder,
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
        self.arch_name = "legacy_mlp"
        if save_path:
            checkpoint = {
                "arch": self.arch_name,
                "input_dim": self.input_dim,
                "hidden_dim": int(self.hidden_dim),
                "num_layers": int(self.num_layers),
                "dropout": float(self.dropout),
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

    def _prepare_input(self, embedding: torch.Tensor) -> torch.Tensor:
        if not torch.is_tensor(embedding):
            embedding = torch.tensor(embedding, dtype=torch.float32)
        vec = embedding.to(self.device).float()
        if vec.dim() == 3:
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
