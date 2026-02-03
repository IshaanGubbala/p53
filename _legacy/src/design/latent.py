from __future__ import annotations

import torch
import torch.nn.functional as F
from transformers import EsmForMaskedLM, EsmTokenizer
from typing import Optional, List, Tuple
from pathlib import Path
import os
from src.core.logging import get_logger

class LatentEmbedder:
    """
    Handles interactions with the ESM-2 Protein Language Model.
    Focuses on encoding sequences into latent space and decoding them back.
    """
    def __init__(self, model_name: str = "facebook/esm2_t6_8M_UR50D", device: Optional[str] = None):
        self.logger = get_logger(__name__)
        self.model_name = model_name
        
        if device is None:
            if torch.backends.mps.is_available():
                self.device = torch.device("mps")
            elif torch.cuda.is_available():
                self.device = torch.device("cuda")
            else:
                self.device = torch.device("cpu")
        else:
            self.device = torch.device(device)
            
        self.logger.info(f"Loading ESM-2 model {model_name} on {self.device}...")
        try:
            self.tokenizer = EsmTokenizer.from_pretrained(model_name)
            self.model = EsmForMaskedLM.from_pretrained(model_name).to(self.device)
            self.model.eval()
            self.logger.info("Model loaded successfully.")
        except Exception as e:
            self.logger.error(f"Failed to load model: {e}")
            raise

    def encode(self, sequence: str) -> torch.Tensor:
        """
        Embeds a protein sequence into the latent space (last hidden state).
        Returns: Tensor of shape (1, L, D) - typically 1280 dim for t33.
        """
        inputs = self.tokenizer(sequence, return_tensors="pt", add_special_tokens=False).to(self.device)
        with torch.no_grad():
            outputs = self.model(**inputs, output_hidden_states=True)
            # Last hidden state is typically used as the "embedding" for per-residue tasks
            # outputs.hidden_states[-1] is the same as outputs.last_hidden_state
            # Shape: (Batch, Seq_Len, Dim)
            embeddings = outputs.hidden_states[-1]
            
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

class ManifoldWalker:
    """
    Performs vector arithmetic and interpolation in the latent space.
    """
    def __init__(self, embedder: LatentEmbedder):
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
