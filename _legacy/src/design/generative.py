from __future__ import annotations

import random
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any

from src.core.logging import get_logger
from src.design.latent import LatentEmbedder, ManifoldWalker
class GenerativeModel(ABC):
    """Abstract base class for generative protein design models."""
    
    def __init__(self, config: dict[str, Any]):
        self.config = config
        self.logger = get_logger(__name__)

    @abstractmethod
    def generate(self, sequence: str, n_samples: int = 10, **kwargs) -> list[str]:
        """
        Generate new sequences based on the input sequence (seed/scaffold).
        
        Args:
            sequence: The starting sequence (e.g. wild type or mutant).
            n_samples: Number of sequences to generate.
            **kwargs: Model specific arguments (masking, diffusion steps, etc).
            
        Returns:
            List of generated amino acid sequences.
        """
        pass

class EvoDiffGenerator(GenerativeModel):
    """
    Interface for EvoDiff (sequence-based diffusion).
    """
    def __init__(self, config: dict[str, Any]):
        super().__init__(config)
        self.model_name = config.get("model_name", "OASis_38M_v2")
        self._load_model()
    
    def _load_model(self):
        self.logger.info(f"Loading EvoDiff model: {self.model_name}")
        # Placeholder for actual model loading
        # from evodiff.pretrained import OA_DM_38M
        # self.model = OA_DM_38M()
        pass

    def generate(self, sequence: str, n_samples: int = 10, **kwargs) -> list[str]:
        self.logger.info(f"Generating {n_samples} samples using EvoDiff (Mock/Stub)")
        # Real implementation would tokenise sequence, apply noise/mask, and reverse diffuse.
        
        # MOCK IMPLEMENTATION:
        # Generate random mutations at variable positions
        results = []
        for _ in range(n_samples):
            # Mutate 1-3 residues randomly
            seq_list = list(sequence)
            n_muts = random.randint(1, 3)
            for _ in range(n_muts):
                idx = random.randint(0, len(seq_list)-1)
                new_aa = random.choice("ACDEFGHIKLMNPQRSTVWY")
                seq_list[idx] = new_aa
            results.append("".join(seq_list))
        
        return results

class RFDiffusionGenerator(GenerativeModel):
    """
    Interface for RFdiffusion (structure-based diffusion).
    """
    def generate(self, sequence: str, n_samples: int = 10, **kwargs) -> list[str]:
        # RFdiffusion typically works on structure (PDB), not just sequence.
        # This wrapper assumes we might pass a PDB path in kwargs or structure.
        return []

class LatentRescueGenerator(GenerativeModel):
    """
    Uses ESM-2 Latent Space Interpolation to find rescue candidates.
    """
    def __init__(self, config: dict[str, Any]):
        super().__init__(config)
        self.model_name = config.get("model_name", "facebook/esm2_t6_8M_UR50D")
        self.embedder = LatentEmbedder(model_name=self.model_name)
        self.walker = ManifoldWalker(self.embedder)
        
    def generate(self, sequence: str, n_samples: int = 10, **kwargs) -> list[str]:
        """
        Generates rescue candidates by interpolating back towards a target (WT) 
        or steering towards high stability.
        
        Args:
            sequence: The starting sequence (Cancer Mutant).
            kwargs: Must contain 'target_sequence' (WT) for interpolation.
        """
        target_seq = kwargs.get("target_sequence")
        if not target_seq:
            self.logger.warning("No target_sequence provided for LatentRescueGenerator. Returning empty.")
            return []
            
        self.logger.info(f"Interpolating from Mutant to WT in latent space. Samples: {n_samples}")
        
        # Interpolate
        # We generate n_samples steps along the path
        candidates = self.walker.interpolate(sequence, target_seq, steps=n_samples)
        
        # Filter out exact start/end duplicates if desired
        unique_candidates = []
        seen = set()
        for cand in candidates:
            if cand not in seen and cand != sequence:
                unique_candidates.append(cand)
                seen.add(cand)
                
        return unique_candidates

def get_generator(name: str, config: dict[str, Any]) -> GenerativeModel:
    if name.lower() == "evodiff":
        return EvoDiffGenerator(config)
    elif name.lower() == "rfdiffusion":
        return RFDiffusionGenerator(config)
    elif name.lower() == "latent":
        return LatentRescueGenerator(config)
    else:
        raise ValueError(f"Unknown generator: {name}")

def extract_mutations(wt_seq: str, gen_seq: str) -> list[str]:
    """Compare generated sequence to WT and return mutation list (e.g. ['A12D'])."""
    if len(wt_seq) != len(gen_seq):
         # If lengths differ, we might need alignment. For now assume fixed length.
         raise ValueError("Sequence lengths differ")
    
    muts = []
    for i, (wt, gen) in enumerate(zip(wt_seq, gen_seq)):
        if wt != gen:
            # 1-indexed position
            muts.append(f"{wt}{i+1}{gen}")
    return muts
