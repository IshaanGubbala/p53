import torch
import numpy as np
from typing import List, Dict
from p53cad.engine.oracle import FunctionalOracle
from p53cad.engine.latent import ManifoldEmbedder

class SaliencyMap:
    def __init__(self, oracle: FunctionalOracle, embedder: ManifoldEmbedder):
        self.oracle = oracle
        self.embedder = embedder
        self.device = oracle.device

    def attribute_residues(self, sequence: str) -> np.ndarray:
        """
        Calculate importance score using Residue Occlusion.
        Impact = |Score(Original) - Score(Occluded)|
        """
        with torch.no_grad():
            # 1. Base Score
            z = self.embedder.encode(sequence) # (1, L, D)
            base_score = self.oracle.predict(z)
            
            L = z.size(1)
            importances = []
            
            # 2. Occlusion Loop
            # We zero out each residue's latent vector one-by-one
            for i in range(L):
                z_occ = z.clone()
                z_occ[:, i, :] = 0.0 # Occlude
                occ_score = self.oracle.predict(z_occ)
                importances.append(abs(base_score - occ_score))
                
            grads = np.array(importances)
            
        # Normalize 0-1
        if grads.max() > grads.min():
            grads = (grads - grads.min()) / (grads.max() - grads.min())
            
        return grads

    def get_top_hotspots(self, sequence: str, top_n: int = 5) -> List[Dict]:
        """Identify residues the model is most sensitive to."""
        scores = self.attribute_residues(sequence)
        top_indices = np.argsort(scores)[::-1][:top_n]
        
        hotspots = []
        for idx in top_indices:
            hotspots.append({
                "pos": int(idx + 1),
                "aa": sequence[idx],
                "importance": float(scores[idx])
            })
        return hotspots
