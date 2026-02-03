from __future__ import annotations

import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
import numpy as np
from pathlib import Path
from typing import List, Tuple
from tqdm import tqdm

from src.core.logging import get_logger
from src.design.latent import LatentEmbedder

class DMSDataset(Dataset):
    def __init__(self, sequences: List[str], scores: List[float], embedder: LatentEmbedder):
        self.sequences = sequences
        self.scores = torch.tensor(scores, dtype=torch.float32)
        self.embedder = embedder
        self.cache = {}
        
    def __len__(self):
        return len(self.sequences)
    
    def __getitem__(self, idx):
        seq = self.sequences[idx]
        score = self.scores[idx]
        
        if seq in self.cache:
            emb = self.cache[seq]
        else:
            # Embed on the fly (slow) or pre-compute
            # For speed in this demo, strict on-the-fly would be slow. 
            # In production, pre-compute all embeddings.
            # Using CPU to avoid clogging GPU memory with single item workflow if not careful
            # unique ID embedding: (1, L, D) -> Mean pool? or use CLS?
            # ESM-2 uses mean of representations usually for seq tasks
            
            # Using embedder.encode returns (1, L, D)
            # We need a fixed size vector for MLP
            
            with torch.no_grad():
                # We need to handle this carefully.
                # Assuming LatentEmbedder can run on CPU for data loading if needed
                # or we just move to device in training loop.
                raw = self.embedder.encode(seq) # (1, L, D)
                emb = raw.mean(dim=1).squeeze(0).cpu() # (D,)
                
            self.cache[seq] = emb
            
        return emb, score

class FunctionalNet(nn.Module):
    """
    Simple MLP Regressor: Embedding (D) -> Function (1)
    """
    def __init__(self, input_dim: int = 320, hidden_dim: int = 128):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.ReLU(),
            nn.Linear(hidden_dim // 2, 1)
        )
        
    def forward(self, x):
        return self.net(x)

class FunctionalOracle:
    def __init__(self, model_path: Path | str = None, input_dim: int = 320): # 320 for esm2_t6
        self.logger = get_logger(__name__)
        self.device = torch.device("mps" if torch.backends.mps.is_available() else "cpu")
        self.model = FunctionalNet(input_dim=input_dim).to(self.device)
        
        if model_path and Path(model_path).exists():
            self.logger.info(f"Loading oracle from {model_path}")
            self.model.load_state_dict(torch.load(model_path, map_location=self.device))
        self.model.eval()

    def train(self, dms_data, embedder: LatentEmbedder, epochs: int = 10, save_path: Path = None):
        self.logger.info(f"Training FunctionalOracle on {len(dms_data)} sequences...")
        self.model.train()
        
        # Prepare Data
        dataset = DMSDataset(dms_data["sequence"].tolist(), dms_data["score"].tolist(), embedder)
        loader = DataLoader(dataset, batch_size=32, shuffle=True)
        
        optimizer = optim.Adam(self.model.parameters(), lr=1e-3)
        criterion = nn.MSELoss()
        
        for epoch in range(epochs):
            total_loss = 0
            for X, y in tqdm(loader, desc=f"Epoch {epoch+1}"):
                X, y = X.to(self.device), y.to(self.device).unsqueeze(1)
                
                optimizer.zero_grad()
                pred = self.model(X)
                loss = criterion(pred, y)
                loss.backward()
                optimizer.step()
                
                total_loss += loss.item()
            
            avg_loss = total_loss / len(loader)
            self.logger.info(f"Epoch {epoch+1} Loss: {avg_loss:.4f}")
            
        self.model.eval()
        if save_path:
            torch.save(self.model.state_dict(), save_path)
            self.logger.info(f"Saved model to {save_path}")

    def predict(self, embedding: torch.Tensor) -> float:
        """
        Predict function score from embedding (1, L, D). 
        """
        with torch.no_grad():
            # Mean pool
            vec = embedding.mean(dim=1) # (1, D)
            return self.model(vec).item()
            
    def predict_batch(self, embeddings: torch.Tensor) -> torch.Tensor:
        """
        Predict function for a batch of embeddings (B, L, D).
        """
        # Embeddings might vary in L if not padded. Assuming pooled input or padded.
        # If input is (B, L, D), we pool:
        vec = embeddings.mean(dim=1)
        return self.model(vec)
