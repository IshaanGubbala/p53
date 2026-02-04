from __future__ import annotations

import torch
import torch.nn.functional as F
from transformers import EsmForMaskedLM, EsmTokenizer
from typing import Optional, List, Tuple
from pathlib import Path
import os
from p53cad.core.logging import get_logger

class ManifoldEmbedder:
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
            # Ensure hidden states are always returned to avoid NoneType errors during navigation
            self.model.config.output_hidden_states = True 
            self.model.eval()
            self.logger.info("Model loaded successfully with output_hidden_states=True.")
        except Exception as e:
            self.logger.error(f"Failed to load model: {e}")
            raise

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

    def latent_forward_ascent(self, embeddings: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        """
        Runs the transformer forward pass starting from soft embeddings.
        Returns: (last_hidden_state, logits, probabilities)
        """
        # Optimized: Explicit hidden state extraction (D=320)
        # Using .esm directly ensures we get the 320-dim latent state for the oracle.
        esm_outputs = self.model.esm(inputs_embeds=embeddings, output_hidden_states=True, return_dict=True)
        h = esm_outputs.last_hidden_state
        
        logits = self.model.lm_head(h)
        probs = torch.softmax(logits[0], dim=-1)
        return h, logits, probs

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
