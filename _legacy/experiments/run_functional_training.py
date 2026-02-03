from __future__ import annotations

import argparse
import sys
import yaml
from pathlib import Path
import pandas as pd

# Add project root
sys.path.append(str(Path.cwd()))

from src.core.logging import get_logger, setup_logging
from src.data.dms_loader import load_dms_data
from src.design.latent import LatentEmbedder
from src.scoring.functional_predictor import FunctionalOracle

def main():
    parser = argparse.ArgumentParser(description="Train Functional Oracle")
    parser.add_argument("--dms_path", type=Path, default="data/raw/p53_DMS_Giacomelli_2018.csv")
    parser.add_argument("--epochs", type=int, default=5)
    parser.add_argument("--model_name", type=str, default="facebook/esm2_t6_8M_UR50D")
    parser.add_argument("--output_dir", type=Path, default="data/models")
    
    args = parser.parse_args()
    setup_logging()
    logger = get_logger(__name__)
    
    # 1. Load Data
    df = load_dms_data(args.dms_path)
    logger.info(f"Loaded {len(df)} records.")
    
    # 2. Fix Sequence
    # Mock data might not have full sequence, just 'mutation' and 'wt'/'alt'.
    # We need the full sequence for ESM.
    # We need a reference WT sequence.
    # Ideally we load WT from fasta and apply mutations.
    
    # Check if df has 'sequence' column
    if "sequence" not in df.columns:
        logger.info("Constructing full sequences from mutations...")
        # Load WT P53
        # Assuming standard location or mocking it to 393 A's if not found (just for fail safe)
        try:
            wt_path = Path("data/interim/uniprot/P04637.fasta")
            if wt_path.exists():
                from experiments.run_design_rescues import _load_sequence
                wt_seq = _load_sequence(wt_path)
            else:
                logger.warning("WT FASTA not found. Using dummy WT for mock training.")
                wt_seq = "M" * 393 # Dummy
                
            # Apply mutations
            sequences = []
            valid_indices = []
            
            for i, row in df.iterrows():
                try:
                    # Parse mutation: e.g. A123V
                    # If mock data, we have 'mutation' col
                    mut = row["mutation"]
                    # If mock data, we also put 'pos', 'alt' in dms_loader
                    
                    # Manual parse if needed
                    import re
                    match = re.match(r"([A-Z])(\d+)([A-Z])", mut)
                    if not match:
                        continue
                        
                    wt_aa, pos, alt_aa = match.groups()
                    pos = int(pos) - 1 # 0-index
                    
                    if pos >= len(wt_seq):
                        continue
                        
                    s_list = list(wt_seq)
                    s_list[pos] = alt_aa
                    seq = "".join(s_list)
                    sequences.append(seq)
                    valid_indices.append(i)
                except Exception:
                    continue
            
            df = df.iloc[valid_indices].copy()
            df["sequence"] = sequences
            logger.info(f"Constructed {len(df)} sequences.")
            
        except Exception as e:
            logger.error(f"Failed to construct sequences: {e}")
            return 1

    # 3. Initialize Embedder
    # t6 has 320 dim
    embedder = LatentEmbedder(model_name=args.model_name)
    
    # 4. Initialize Oracle
    oracle = FunctionalOracle(input_dim=320) # ensure matches model
    
    # 5. Train
    args.output_dir.mkdir(parents=True, exist_ok=True)
    save_path = args.output_dir / "functional_oracle.pt"
    
    try:
        oracle.train(df, embedder, epochs=args.epochs, save_path=save_path)
    except KeyboardInterrupt:
        logger.warning("Training interrupted. Saving current state...")
        import torch
        torch.save(oracle.model.state_dict(), save_path)
        
    logger.info("Training complete.")
    return 0

if __name__ == "__main__":
    main()
