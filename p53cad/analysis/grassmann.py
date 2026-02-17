from __future__ import annotations

import torch
import torch.nn.functional as F
from p53cad.engine.latent import ManifoldEmbedder
from p53cad.core.logging import get_logger
import numpy as np

class GrassmannMetric:
    """
    Analyzes the 'Functional Geometry' of proteins by treating them as subspaces on a Grassmann manifold.
    
    The core idea:
    A protein's function is encoded in how its residues attend to each other.
    We extract the 'Key' or 'Value' matrices from the Multi-Head Attention layers.
    These define a subspace in R^D.
    We compare the subspace of a Mutant vs Wild-Types using Principal Angles.
    """
    def __init__(self, embedder: ManifoldEmbedder):
        self.embedder = embedder
        self.logger = get_logger(__name__)

    def extract_subspace(self, sequence: str, layer: int = -1) -> torch.Tensor:
        """
        Extracts the subspace representation (U) of the sequence from a specific layer.
        Ideally, we want the internal attention matrices, but extracting those from 
        HuggingFace ESM model requires hooks or output_attentions=True.
        
        Simplified Proxy:
        We use the Hidden States (L, D) and perform SVD.
        The top-k singular vectors define the 'Activity Subspace' of the protein.
        """
        # (1, L, D)
        z = self.embedder.encode(sequence)
        
        # Squeeze to (L, D)
        matrix = z.squeeze(0)
        
        # SVD on CPU to avoid MPS 'aten::_linalg_svd.U' NotImplementedError
        # The matrix is (L, D), typically (393, 320) or similar.
        # This is small enough that CPU compute is negligible.
        cpu_matrix = matrix.detach().cpu()
        U, S, Vh = torch.linalg.svd(cpu_matrix, full_matrices=False)
        # Vh is (min(L,D), D). The rows are the singular vectors.
        # Move back to original device if needed, or keep on CPU for now.
        # Grassmann distance logic works on whatever device the subspaces are on.
        
        # Let's take top 10 principal components
        k = 10
        if Vh.shape[0] < k:
            k = Vh.shape[0]
            
        subspace = Vh[:k, :].to(matrix.device) # (k, D)
        # Orthonormalize just in case (SVD gives orthonormal usually)
        return subspace

    def grassmann_distance(self, seq_a: str, seq_b: str) -> float:
        """
        Computes the Grassmann Distance (Geodesic distance on Gr(k, D)) between two sequences.
        d(A, B) = sqrt(sum(theta_i^2))
        where theta_i are the principal angles.
        """
        sub_a = self.extract_subspace(seq_a)
        sub_b = self.extract_subspace(seq_b)
        
        # Compute Principal Angles
        # Bjorck & Golub algorithm or just simple SVD of product
        # If A and B are orthonormal bases (k, D):
        # M = A @ B.T  --> (k, k)
        # Cosines are singular values of M.
        
        # A and B (sub_a, sub_b) are (k, D)
        # We need to transpose to match A @ B.T if they were (D, k).
        # here they are (k, D). So A @ B.T is (k, k).
        
        # Ensure on same device
        if sub_a.device != sub_b.device:
            sub_b = sub_b.to(sub_a.device)
            
        interaction = torch.matmul(sub_a, sub_b.T)
        
        # Singular values of interaction matrix are the cosines of principal angles
        # range [0, 1]
        # Move to CPU for SVD because 'aten::_linalg_svd.U' is not implemented on MPS
        s_vals = torch.linalg.svdvals(interaction.cpu())
        
        # Clamp for numerical safety
        s_vals = torch.clamp(s_vals, -1.0, 1.0)
        
        # Angles
        thetas = torch.acos(s_vals)
        
        # Geodesic distance (L2 norm of angles)
        dist = torch.norm(thetas).item()
        
        return dist
