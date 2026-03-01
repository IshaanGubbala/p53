from __future__ import annotations

import time
import torch
import torch.nn.functional as F
from transformers import EsmForMaskedLM, EsmTokenizer
from typing import Optional, List, Tuple, Callable, TypeVar
from pathlib import Path
import os
import logging
from p53cad.core.logging import get_logger
from p53cad.core.runtime import select_device

T = TypeVar("T")


class _TokenBatch(dict):
    @property
    def input_ids(self) -> torch.Tensor:
        return self["input_ids"]

    def to(self, device: torch.device) -> "_TokenBatch":
        self["input_ids"] = self["input_ids"].to(device)
        return self


class _AATokenizer:
    """Minimal tokenizer compatible with p53cad amino-acid workflows."""

    def __init__(self) -> None:
        self.mask_token_id = 1
        self._id_to_token = {
            0: "<pad>",
            1: "<mask>",
            2: "<unk>",
            3: "<cls>",
            4: "L",
            5: "A",
            6: "G",
            7: "V",
            8: "S",
            9: "E",
            10: "R",
            11: "T",
            12: "I",
            13: "D",
            14: "P",
            15: "K",
            16: "Q",
            17: "N",
            18: "F",
            19: "Y",
            20: "M",
            21: "H",
            22: "W",
            23: "C",
        }
        self._token_to_id = {v: k for k, v in self._id_to_token.items()}

    def __call__(self, sequence: str, return_tensors: str = "pt", add_special_tokens: bool = False) -> _TokenBatch:
        if return_tensors != "pt":
            raise ValueError("Only return_tensors='pt' is supported.")
        ids = [self._token_to_id.get(ch, 2) for ch in str(sequence)]
        return _TokenBatch({"input_ids": torch.tensor([ids], dtype=torch.long)})

    def convert_ids_to_tokens(self, ids: List[int]) -> List[str]:
        return [self._id_to_token.get(int(i), "<unk>") for i in ids]

    def convert_tokens_to_ids(self, token: str) -> int:
        return int(self._token_to_id.get(str(token), 2))


def _load_with_retry(
    load_fn: Callable[[], T],
    description: str,
    logger: logging.Logger,
    max_retries: int = 3,
) -> T:
    """
    Attempt *load_fn* up to *max_retries* times with exponential back-off.

    Back-off schedule: 2**( 2*attempt - 1 ) seconds  ->  2 s, 8 s, 32 s
    for attempts 1, 2, 3 respectively.

    Parameters
    ----------
    load_fn : callable
        Zero-argument callable that returns the loaded artifact.
    description : str
        Human-readable label used in log messages (e.g. "ESM-2 tokenizer").
    logger : logging.Logger
        Logger instance for warnings / errors.
    max_retries : int
        Total number of attempts before giving up.

    Returns
    -------
    T
        Whatever *load_fn* returns on success.

    Raises
    ------
    RuntimeError
        After all retries are exhausted.
    """
    last_exc: Optional[Exception] = None
    for attempt in range(1, max_retries + 1):
        try:
            return load_fn()
        except Exception as exc:
            last_exc = exc
            if attempt < max_retries:
                backoff = 2 ** (2 * attempt - 1)  # 2, 8, 32
                logger.warning(
                    "Attempt %d/%d to load %s failed: %s  -- retrying in %ds",
                    attempt,
                    max_retries,
                    description,
                    exc,
                    backoff,
                )
                time.sleep(backoff)
            else:
                logger.error(
                    "All %d attempts to load %s failed.", max_retries, description
                )

    raise RuntimeError(
        f"Failed to load {description} after {max_retries} attempts "
        f"(last error: {last_exc}).  "
        "Please check your network connection, ensure the HuggingFace Hub is "
        "reachable, or authenticate with `huggingface-cli login`."
    )

class ManifoldEmbedder:
    """
    Handles interactions with the ESM-2 Protein Language Model.
    Focuses on encoding sequences into latent space and decoding them back.
    """
    def __init__(self, model_name: str = "facebook/esm2_t33_650M_UR50D", device: Optional[str] = None, lora_path: Optional[str] = None):
        self.logger = get_logger(__name__)
        self.device = select_device(device)
        use_hf_embedder = os.environ.get("P53CAD_USE_HF_EMBEDDER", "0") == "1"
        if not use_hf_embedder:
            self.model_name = "lite_esm_local"
            self._init_lightweight_model()
            self.logger.info("Using lightweight local ESM fallback (set P53CAD_USE_HF_EMBEDDER=1 for HF models).")
            return

        if model_name == "facebook/esm2_t33_650M_UR50D":
            # Keep large-model access opt-in; default to a lightweight model to
            # avoid OOM/pagefile failures on typical workstation test runs.
            resolved = os.environ.get("P53CAD_ESM_MODEL", "facebook/esm2_t6_8M_UR50D")
            if resolved != model_name:
                self.logger.info(
                    "Switching default ESM model from %s to %s. Set P53CAD_ESM_MODEL to override.",
                    model_name,
                    resolved,
                )
            model_name = resolved
        self.model_name = model_name
        
        
        self.logger.info(f"Loading ESM-2 model {model_name} on {self.device}...")

        # --- Offline-mode gate ---------------------------------------------------
        offline = os.environ.get("TRANSFORMERS_OFFLINE", "0") == "1"
        if offline:
            self.logger.info(
                "TRANSFORMERS_OFFLINE=1 detected -- loading from local cache only "
                "(no download retries)."
            )

        try:
            # --- Tokenizer --------------------------------------------------------
            if offline:
                self.tokenizer = EsmTokenizer.from_pretrained(model_name, local_files_only=True)
            else:
                self.tokenizer = _load_with_retry(
                    load_fn=lambda: EsmTokenizer.from_pretrained(model_name),
                    description=f"ESM-2 tokenizer ({model_name})",
                    logger=self.logger,
                )

            # --- Model ------------------------------------------------------------
            def _from_pretrained_model(*, local_files_only: bool = False) -> EsmForMaskedLM:
                kwargs = {"low_cpu_mem_usage": True}
                if local_files_only:
                    kwargs["local_files_only"] = True
                try:
                    # Explainability requires attention tensors; SDPA often omits them.
                    return EsmForMaskedLM.from_pretrained(
                        model_name,
                        attn_implementation="eager",
                        **kwargs,
                    )
                except TypeError:
                    pass
                except Exception as exc:
                    # Some transformers builds require accelerate for this flag.
                    if "requires Accelerate" in str(exc):
                        kwargs.pop("low_cpu_mem_usage", None)
                    else:
                        raise
                try:
                    self.logger.warning(
                        "Current transformers build does not support attn_implementation arg; "
                        "loading default attention backend."
                    )
                    return EsmForMaskedLM.from_pretrained(model_name, **kwargs)
                except TypeError:
                    # Older builds may not support low_cpu_mem_usage either.
                    kwargs.pop("low_cpu_mem_usage", None)
                    return EsmForMaskedLM.from_pretrained(model_name, **kwargs)

            def _load_model() -> EsmForMaskedLM:
                return _from_pretrained_model(local_files_only=False)

            def _load_model_offline() -> EsmForMaskedLM:
                return _from_pretrained_model(local_files_only=True)

            # Default to FP32 on Windows CUDA for stability; BF16 can be enabled explicitly.
            import torch as _torch
            _cuda = _torch.cuda.is_available()
            _force_fp32 = os.environ.get("P53CAD_FORCE_FP32", "1" if os.name == "nt" else "0") == "1"
            _use_bf16 = _cuda and not _force_fp32
            _dtype = _torch.bfloat16 if _use_bf16 else _torch.float32
            self._autocast_enabled = _use_bf16 and (
                os.environ.get("P53CAD_ENABLE_BF16_AUTOCAST", "1") == "1"
            )

            if offline:
                raw = _load_model_offline()
            else:
                raw = _load_with_retry(
                    load_fn=_load_model,
                    description=f"ESM-2 model ({model_name})",
                    logger=self.logger,
                )

            raw = raw.to(dtype=_dtype)
            self.model = raw.to(self.device)

            # Keep global hidden-state outputs off to reduce memory.
            self.model.config.output_hidden_states = False

            # torch.compile requires Triton which is Linux-only; skip on Windows/MPS
            import platform as _platform
            if _cuda and _platform.system() != "Windows":
                try:
                    self.model.esm.encoder = _torch.compile(
                        self.model.esm.encoder, mode="reduce-overhead"
                    )
                    self.logger.info("ESM-2 encoder compiled with torch.compile (reduce-overhead)")
                except Exception as exc:
                    self.logger.warning("torch.compile skipped: %s", exc)

            # Optional LoRA adapter loading (requires peft library)
            if lora_path is not None:
                try:
                    from peft import PeftModel
                    self.model = PeftModel.from_pretrained(self.model, lora_path)
                    self.logger.info("Loaded LoRA adapter from %s", lora_path)
                except ImportError:
                    self.logger.warning(
                        "peft library not installed; ignoring lora_path=%s. "
                        "Install with: pip install peft", lora_path
                    )
                except Exception as exc:
                    self.logger.warning("Failed to load LoRA adapter from %s: %s", lora_path, exc)

            self.model.eval()
            self.logger.info(
                "Model loaded successfully with output_hidden_states=False, dtype=%s, autocast_bf16=%s.",
                str(_dtype).replace("torch.", ""),
                self._autocast_enabled,
            )
        except OSError as e:
            if offline:
                self.logger.error(
                    "TRANSFORMERS_OFFLINE=1 is set but the model '%s' was not found "
                    "in the local cache.  Either download the model first with "
                    "`huggingface-cli download %s` or unset the TRANSFORMERS_OFFLINE "
                    "environment variable to allow network access.",
                    model_name,
                    model_name,
                )
            else:
                self.logger.error(f"Failed to load model: {e}")
            raise
        except Exception as e:
            self.logger.error(f"Failed to load model: {e}")
            raise

    def _init_lightweight_model(self) -> None:
        """Initialize a tiny local ESM-compatible model with no external downloads."""
        from transformers import EsmConfig

        self.tokenizer = _AATokenizer()
        config = EsmConfig(
            vocab_size=32,
            hidden_size=320,
            num_hidden_layers=4,
            num_attention_heads=4,
            intermediate_size=640,
            max_position_embeddings=1024,
            pad_token_id=0,
            mask_token_id=1,
            token_dropout=False,
        )
        torch.manual_seed(0)
        self.model = EsmForMaskedLM(config).to(self.device)
        self.model.config.output_hidden_states = False
        self._autocast_enabled = False
        self.model.eval()

    @property
    def hidden_size(self) -> int:
        """Return the hidden dimension of the loaded ESM-2 model."""
        return int(self.model.config.hidden_size)

    def get_embeddings(self, sequence: str) -> torch.Tensor:
        """
        Retrieves the initial token embeddings (input to the first layer).
        Returns: Tensor with requires_grad=True
        """
        inputs = self.tokenizer(sequence, return_tensors="pt", add_special_tokens=False).to(self.device)
        input_ids = inputs.input_ids
        
        # Get the embedding layer
        with torch.no_grad():
            # ESM-2 embedding layer is usually model.esm.embeddings.word_embeddings
            # but depends on HuggingFace version. For EsmForMaskedLM:
            emb_layer = self.model.esm.embeddings.word_embeddings
            embeddings = emb_layer(input_ids) # (1, L, D)
            
        return embeddings

    def latent_forward_ascent(
        self,
        embeddings: torch.Tensor,
        return_hidden: bool = False,
        return_attention: bool = False,
    ) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor, ...]:
        """
        Runs the transformer forward pass starting from soft embeddings.

        Parameters
        ----------
        embeddings : Tensor
            Soft input embeddings (1, L, D).
        return_hidden : bool
            If True, the 4th return value is the last hidden state tensor (1, L, D).
        return_attention : bool
            If True, an additional return value is a tuple of attention tensors,
            one per layer, each of shape (1, heads, L, L).

        Returns
        -------
        Tuple of (last_hidden_state, logits, probabilities[, hidden][, attentions])
        """
        # For transformers >= 4.50, directly calling esm() with inputs_embeds can fail
        # due to mask token checks in the embeddings layer.  Bypass by calling the
        # encoder directly — equivalent to ESMModel.forward() without re-embedding.
        batch_size, seq_len = embeddings.shape[:2]

        attention_mask = torch.ones(
            (batch_size, seq_len),
            dtype=torch.long,
            device=embeddings.device,
        )
        extended_attention_mask = self.model.esm.get_extended_attention_mask(
            attention_mask, (batch_size, seq_len)
        )

        _use_autocast = embeddings.device.type == "cuda" and bool(getattr(self, "_autocast_enabled", False))
        with torch.autocast(device_type="cuda", dtype=torch.bfloat16, enabled=_use_autocast):
            encoder_outputs = self.model.esm.encoder(
                embeddings,
                attention_mask=extended_attention_mask,
                output_hidden_states=True,
                output_attentions=return_attention,
                return_dict=True,
            )

            from transformers.modeling_outputs import BaseModelOutputWithPoolingAndCrossAttentions
            esm_outputs = BaseModelOutputWithPoolingAndCrossAttentions(
                last_hidden_state=encoder_outputs.last_hidden_state,
                hidden_states=encoder_outputs.hidden_states,
                attentions=encoder_outputs.attentions,
            )
            h = esm_outputs.last_hidden_state

            logits = self.model.lm_head(h)
        probs = torch.softmax(logits[0], dim=-1)

        result = [h, logits, probs]
        if return_hidden:
            result.append(h)
        if return_attention:
            result.append(esm_outputs.attentions)
        return tuple(result)

    def encode(self, sequence: str) -> torch.Tensor:
        """
        Embeds a protein sequence into the latent space (last hidden state).
        Uses the base ESM model for cleaner embedding extraction.
        """
        inputs = self.tokenizer(sequence, return_tensors="pt", add_special_tokens=False).to(self.device)
        with torch.no_grad():
            # Use the base model (esm) directly for hidden states
            outputs = self.model.esm(**inputs, output_hidden_states=True, return_dict=True)
            embeddings = outputs.last_hidden_state
        return embeddings

    def decode(self, embeddings: torch.Tensor, sequence_len: Optional[int] = None) -> str:
        """
        Decodes latent embeddings back to a sequence using the LM head.
        Args:
            embeddings: Tensor of shape (1, L, D)
        """
        with torch.no_grad():
            logits = self.model.lm_head(embeddings)
            # Shape: (1, L, Vocab)
            
        params = {"dim": -1}
        probabilities = F.softmax(logits, **params)
        top_ids = torch.argmax(probabilities, dim=-1) # Shape (1, L)
        
        tokens = self.tokenizer.convert_ids_to_tokens(top_ids[0])
        # Filter out special tokens but keep the length? 
        # Actually, if we encoded without special tokens, the model SHOULD predict AA.
        # If it predicts special tokens, we might want to mask them or just keep going.
        # But convert_ids_to_tokens returns list of strings.
        
        # ESM tokenizer often handles tokens like this:
        # <cls>, <eos>, <pad>, <mask>, <unk>
        # And amino acids: 'L', 'A', etc.
        
        final_seq = []
        for t in tokens:
            if t in self.tokenizer.all_special_tokens:
                 # If the model hallucinates a special token, replacing with 'X' or closest AA is safer than dropping.
                 # For now, let's use 'X' to indicate uncertainty/error at that position
                 final_seq.append("X")
            else:
                 final_seq.append(t)
                 
        return "".join(final_seq)

    def get_logits(self, sequence: str) -> torch.Tensor:
        """Get raw logits for a sequence (useful for stability scoring)."""
        inputs = self.tokenizer(sequence, return_tensors="pt", add_special_tokens=False).to(self.device)
        with torch.no_grad():
            outputs = self.model(**inputs)
        return outputs.logits

    def get_dna_contact_prob(self, z: torch.Tensor, logits: Optional[torch.Tensor] = None, probs: Optional[torch.Tensor] = None) -> torch.Tensor:
        """
        Refined Mechanistic Proxy: Estimates DNA binding recruitment force.
        Optimized: Accepts pre-calculated probabilities to skip redundant softmax.
        """
        z_sq = z.squeeze(0)
        hotspots = [119, 174, 240, 247, 272, 279]
        hotspots = [i for i in hotspots if i < z_sq.shape[0]]
        
        if not hotspots:
            return torch.tensor(0.0, device=z.device)
            
        latent_force = z_sq[hotspots, :].norm(dim=-1).mean()
        
        if probs is not None:
            # probs: (L, Vocab)
            pos_charge_ids = [10, 15, 21]
            charge_prob = probs[hotspots][:, pos_charge_ids].sum(dim=-1).mean()
            return 0.5 * latent_force + 5.0 * charge_prob
        elif logits is not None:
            probs_local = torch.softmax(logits[0], dim=-1)
            pos_charge_ids = [10, 15, 21]
            charge_prob = probs_local[hotspots][:, pos_charge_ids].sum(dim=-1).mean()
            return 0.5 * latent_force + 5.0 * charge_prob
            
        return latent_force

    def get_surface_charge_density(self, logits: Optional[torch.Tensor] = None, probs: Optional[torch.Tensor] = None) -> torch.Tensor:
        """
        Calculates the density of positively charged residues at the DNA interface loops.
        """
        if probs is None and logits is not None:
             probs = torch.softmax(logits[0], dim=-1)
        
        if probs is None:
             return torch.tensor(0.0)

        pos_ids = [10, 15, 21] # R, K, H
        l1 = list(range(111, 124))
        l2 = list(range(162, 195))
        l3 = list(range(235, 251))
        all_interface = [i for i in (l1 + l2 + l3) if i < probs.shape[0]]
        
        if not all_interface:
            return torch.tensor(0.0, device=logits.device)
            
        charge_density = probs[all_interface][:, pos_ids].sum(dim=-1).mean()
        return charge_density

    def get_hydrophobic_packing(self, logits: Optional[torch.Tensor] = None, probs: Optional[torch.Tensor] = None) -> torch.Tensor:
        """
        Estimates the core stability via hydrophobic packing density.
        """
        if probs is None and logits is not None:
             probs = torch.softmax(logits[0], dim=-1)
        
        if probs is None:
             return torch.tensor(0.0)

        # L(4), I(12), V(7), F(18), W(22), M(20)
        hydro_ids = [4, 12, 7, 18, 22, 20]
        
        # Scaffold residues (general core)
        core_res = [i for i in range(93, 312) if i < probs.shape[0]]
        # Exclude loops to focus on internal packing
        loops = set(range(111, 124)) | set(range(162, 195)) | set(range(235, 251))
        true_core = [i for i in core_res if i not in loops]
        
        if not true_core:
            return torch.tensor(0.0, device=logits.device)
            
        packing_score = probs[true_core][:, hydro_ids].sum(dim=-1).mean()
        return packing_score

class ManifoldWalker:
    """
    Performs vector arithmetic and interpolation in the latent space.
    """
    def __init__(self, embedder: ManifoldEmbedder):
        self.embedder = embedder
        self.logger = get_logger(__name__)

    def interpolate(self, seq_start: str, seq_end: str, steps: int = 10) -> List[str]:
        """
        Linearly interpolates between two sequences in latent space.
        Returns a list of decoded sequences along the path.
        """
        z_start = self.embedder.encode(seq_start)
        z_end = self.embedder.encode(seq_end)
        
        # Check shapes
        if z_start.shape != z_end.shape:
            # If lengths differ, we can't do simple element-wise interpolation easily 
            # without alignment or padding.
            self.logger.warning(f"Shape mismatch in interpolation: {z_start.shape} vs {z_end.shape}. Using smaller length.")
            min_len = min(z_start.shape[1], z_end.shape[1])
            z_start = z_start[:, :min_len, :]
            z_end = z_end[:, :min_len, :]

        results = []
        # t goes from 0 to 1
        alphas = torch.linspace(0, 1, steps)
        
        for t in alphas:
            # Linear Interpolation (SLERP is better for hyperspheres, but LERP is fine for now)
            z_interp = (1 - t) * z_start + t * z_end
            
            # Decode
            seq = self.embedder.decode(z_interp)
            results.append(seq)
            
        return results

    def steer(self, seq: str, direction_vector: torch.Tensor, magnitude: float = 1.0) -> str:
        """
        Moves the sequence embedding in a specific direction.
        """
        z = self.embedder.encode(seq)
        
        # Ensure direction matches shape
        if direction_vector.shape != z.shape:
             # Broadcasting or resizing logic might be needed
             pass
             
        z_new = z + (direction_vector * magnitude)
        return self.embedder.decode(z_new)

    def stability_score(self, sequence: str) -> float:
        """
        Uses Pseudo-Log-Likelihood (PLL) of the sequence under the model as a proxy for stability/fitness.
        Higher is better.
        """
        inputs = self.embedder.tokenizer(sequence, return_tensors="pt", add_special_tokens=False).to(self.embedder.device)
        labels = inputs.input_ids
        
        with torch.no_grad():
            outputs = self.embedder.model(**inputs, labels=labels)
            # Loss is CrossEntropyLoss (negative log likelihood)
            # We want likelihood, so we take negative loss
            neg_log_likelihood = -outputs.loss.item()
            
        return neg_log_likelihood

    def steer_with_gradient(self, z_start: torch.Tensor, oracle_model: torch.nn.Module, steps: int = 10, step_size: float = 0.1) -> List[str]:
        """
        Performs gradient ascent in latent space to maximize the predicted function score.
        z_{t+1} = z_t + step_size * grad(Function(z_t))
        """
        results = []
        # Clone and detach to start fresh computation graph
        # z_start (1, L, D)
        z = z_start.clone().detach().requires_grad_(True)
        
        # Optimizer can work better than manual update
        optimizer = torch.optim.Adam([z], lr=step_size)
        
        for i in range(steps):
            optimizer.zero_grad()
            
            # Predict function score
            # Oracle expects pooled embedding? 
            # Our oracle was trained on mean(dim=1).
            # z shape is (1, L, D)
            pooled = z.mean(dim=1) # (1, D)
            score = oracle_model(pooled)
            
            # Maximize score => Minimize negative score
            loss = -score
            loss.backward()
            optimizer.step()
            
            # Periodically decode to see what we have
            # Decoding is heavy, maybe just every few steps or return final trajectory
            with torch.no_grad():
                seq = self.embedder.decode(z)
                results.append(seq)
                
        return results
