"""LoRA fine-tuning of ESM-2 on p53 DMS data.

Provides domain-specific adaptation of the pretrained ESM-2 model using
Low-Rank Adaptation (LoRA) via the ``peft`` library.  Two training modes:

1. **Contrastive training**: Triplet loss where anchor=WT, positive=functional
   variants, negative=LoF variants.  Pulls functional rescue embeddings
   together and pushes LoF embeddings apart.

2. **Masked language model (MLM)**: Domain adaptation on p53 family sequences
   to improve the model's understanding of p53-specific sequence constraints.

The ``peft`` library is optional — this module raises ``ImportError`` at
construction time if it is not installed.

Usage::

    from p53cad.engine.finetune import ESM2FineTuner

    finetuner = ESM2FineTuner("facebook/esm2_t33_650M_UR50D")
    finetuner.train_contrastive(functional_seqs, lof_seqs, epochs=5)
    finetuner.save_adapter("data/models/esm2_lora_p53")
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import List, Optional

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader, Dataset

from p53cad.core.logging import get_logger

logger = get_logger(__name__)


class _TripletDataset(Dataset):
    """Dataset for contrastive triplet training."""

    def __init__(self, anchor_ids, positive_ids, negative_ids):
        self.anchor_ids = anchor_ids
        self.positive_ids = positive_ids
        self.negative_ids = negative_ids

    def __len__(self):
        return len(self.anchor_ids)

    def __getitem__(self, idx):
        return self.anchor_ids[idx], self.positive_ids[idx], self.negative_ids[idx]


class _MLMDataset(Dataset):
    """Dataset for masked language model training."""

    def __init__(self, token_ids_list, mask_prob=0.15, mask_token_id=32):
        self.token_ids_list = token_ids_list
        self.mask_prob = mask_prob
        self.mask_token_id = mask_token_id

    def __len__(self):
        return len(self.token_ids_list)

    def __getitem__(self, idx):
        ids = self.token_ids_list[idx].clone()
        labels = ids.clone()

        # Randomly mask positions
        mask = torch.rand(ids.shape) < self.mask_prob
        ids[mask] = self.mask_token_id
        labels[~mask] = -100  # ignore non-masked positions in loss

        return ids, labels


class ESM2FineTuner:
    """LoRA fine-tuning of ESM-2 on p53 functional data.

    Parameters
    ----------
    model_name : str
        HuggingFace model identifier.
    lora_r : int
        LoRA rank (default 8).
    lora_alpha : int
        LoRA scaling factor (default 16).
    target_modules : list of str
        Attention weight matrices to adapt.
    device : str or None
        Compute device.  Auto-detected if None.
    """

    def __init__(
        self,
        model_name: str = "facebook/esm2_t33_650M_UR50D",
        lora_r: int = 8,
        lora_alpha: int = 16,
        target_modules: Optional[List[str]] = None,
        device: Optional[str] = None,
    ):
        try:
            from peft import LoraConfig, get_peft_model, TaskType
        except ImportError as exc:
            raise ImportError(
                "ESM2FineTuner requires the 'peft' library. "
                "Install with: pip install peft"
            ) from exc

        from transformers import EsmForMaskedLM, EsmTokenizer

        self.model_name = model_name
        self.device = torch.device(
            device
            or ("mps" if torch.backends.mps.is_available() else "cpu")
        )

        logger.info("Loading ESM-2 for LoRA fine-tuning: %s on %s", model_name, self.device)
        self.tokenizer = EsmTokenizer.from_pretrained(model_name)
        self.base_model = EsmForMaskedLM.from_pretrained(model_name)

        if target_modules is None:
            target_modules = ["query", "key", "value"]

        lora_config = LoraConfig(
            r=lora_r,
            lora_alpha=lora_alpha,
            target_modules=target_modules,
            lora_dropout=0.1,
            bias="none",
            task_type=TaskType.FEATURE_EXTRACTION,
        )

        self.model = get_peft_model(self.base_model, lora_config)
        self.model.to(self.device)

        trainable = sum(p.numel() for p in self.model.parameters() if p.requires_grad)
        total = sum(p.numel() for p in self.model.parameters())
        logger.info(
            "LoRA applied: %d trainable / %d total parameters (%.2f%%)",
            trainable, total, 100.0 * trainable / max(total, 1),
        )

    def _encode_sequences(self, sequences: List[str]) -> torch.Tensor:
        """Tokenize sequences and return input_ids tensor."""
        encoded = self.tokenizer(
            sequences, return_tensors="pt", padding=True,
            truncation=True, add_special_tokens=False,
        )
        return encoded.input_ids

    def _get_embeddings(self, input_ids: torch.Tensor) -> torch.Tensor:
        """Get mean-pooled embeddings from model."""
        input_ids = input_ids.to(self.device)
        outputs = self.model.base_model.esm(input_ids=input_ids, output_hidden_states=True, return_dict=True)
        hidden = outputs.last_hidden_state  # (batch, L, D)
        return hidden.mean(dim=1)  # (batch, D)

    def train_contrastive(
        self,
        functional_seqs: List[str],
        lof_seqs: List[str],
        wt_sequence: Optional[str] = None,
        epochs: int = 5,
        batch_size: int = 4,
        lr: float = 1e-4,
        margin: float = 1.0,
    ) -> None:
        """Train with triplet loss: anchor=WT, positive=functional, negative=LoF.

        Parameters
        ----------
        functional_seqs : list of str
            Protein sequences of functional variants (positive examples).
        lof_seqs : list of str
            Protein sequences of loss-of-function variants (negative examples).
        wt_sequence : str or None
            Wild-type sequence used as anchor.  If None, uses p53 WT.
        epochs : int
            Number of training epochs.
        batch_size : int
            Training batch size.
        lr : float
            Learning rate.
        margin : float
            Triplet loss margin.
        """
        if wt_sequence is None:
            from p53cad.data.dms import P53_WT
            wt_sequence = P53_WT

        n_pairs = min(len(functional_seqs), len(lof_seqs))
        if n_pairs == 0:
            logger.warning("No training pairs available for contrastive training")
            return

        logger.info(
            "Starting contrastive LoRA training: %d pairs, %d epochs",
            n_pairs, epochs,
        )

        # Tokenize
        anchor_ids = self._encode_sequences([wt_sequence] * n_pairs)
        positive_ids = self._encode_sequences(functional_seqs[:n_pairs])
        negative_ids = self._encode_sequences(lof_seqs[:n_pairs])

        dataset = _TripletDataset(anchor_ids, positive_ids, negative_ids)
        loader = DataLoader(dataset, batch_size=batch_size, shuffle=True)

        triplet_loss_fn = nn.TripletMarginLoss(margin=margin)
        optimizer = torch.optim.AdamW(
            (p for p in self.model.parameters() if p.requires_grad),
            lr=lr,
        )

        self.model.train()
        for epoch in range(epochs):
            total_loss = 0.0
            n_batches = 0
            for anchor, positive, negative in loader:
                optimizer.zero_grad()
                emb_a = self._get_embeddings(anchor)
                emb_p = self._get_embeddings(positive)
                emb_n = self._get_embeddings(negative)

                loss = triplet_loss_fn(emb_a, emb_p, emb_n)
                loss.backward()
                optimizer.step()

                total_loss += float(loss.item())
                n_batches += 1

            avg_loss = total_loss / max(n_batches, 1)
            logger.info("Contrastive epoch %d/%d: loss=%.4f", epoch + 1, epochs, avg_loss)

        self.model.eval()
        logger.info("Contrastive training complete")

    def train_mlm(
        self,
        sequences: List[str],
        epochs: int = 3,
        batch_size: int = 4,
        lr: float = 5e-5,
        mask_prob: float = 0.15,
    ) -> None:
        """Masked language model training on p53 family sequences.

        Parameters
        ----------
        sequences : list of str
            p53 family protein sequences for domain adaptation.
        epochs : int
            Number of training epochs.
        batch_size : int
            Training batch size.
        lr : float
            Learning rate.
        mask_prob : float
            Fraction of tokens to mask per sequence.
        """
        if not sequences:
            logger.warning("No sequences for MLM training")
            return

        logger.info("Starting MLM LoRA training: %d sequences, %d epochs", len(sequences), epochs)

        token_ids_list = []
        for seq in sequences:
            ids = self.tokenizer(seq, return_tensors="pt", add_special_tokens=False).input_ids.squeeze(0)
            token_ids_list.append(ids)

        dataset = _MLMDataset(
            token_ids_list,
            mask_prob=mask_prob,
            mask_token_id=self.tokenizer.mask_token_id,
        )
        loader = DataLoader(dataset, batch_size=batch_size, shuffle=True)

        optimizer = torch.optim.AdamW(
            (p for p in self.model.parameters() if p.requires_grad),
            lr=lr,
        )

        self.model.train()
        for epoch in range(epochs):
            total_loss = 0.0
            n_batches = 0
            for input_ids, labels in loader:
                input_ids = input_ids.to(self.device)
                labels = labels.to(self.device)

                optimizer.zero_grad()
                outputs = self.model(input_ids=input_ids, labels=labels, return_dict=True)
                loss = outputs.loss
                loss.backward()
                optimizer.step()

                total_loss += float(loss.item())
                n_batches += 1

            avg_loss = total_loss / max(n_batches, 1)
            logger.info("MLM epoch %d/%d: loss=%.4f", epoch + 1, epochs, avg_loss)

        self.model.eval()
        logger.info("MLM training complete")

    def save_adapter(self, path: str | Path) -> Path:
        """Save LoRA adapter weights (not full model)."""
        path = Path(path)
        path.mkdir(parents=True, exist_ok=True)
        self.model.save_pretrained(str(path))
        logger.info("LoRA adapter saved to %s", path)
        return path

    def merge_and_export(self, path: str | Path) -> Path:
        """Merge LoRA weights into base model and save full checkpoint."""
        path = Path(path)
        path.mkdir(parents=True, exist_ok=True)
        merged = self.model.merge_and_unload()
        merged.save_pretrained(str(path))
        self.tokenizer.save_pretrained(str(path))
        logger.info("Merged model saved to %s", path)
        return path
