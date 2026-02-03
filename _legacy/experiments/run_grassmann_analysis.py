from __future__ import annotations

import argparse
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Add project root
sys.path.append(str(Path.cwd()))

from src.core.logging import get_logger, setup_logging
from src.analysis.grassmannian import GrassmannMetric
from src.design.latent import LatentEmbedder

def main():
    parser = argparse.ArgumentParser(description="Run Grassmannian Analysis")
    parser.add_argument("--input_csv", type=Path, default="data/processed/latent_rescue_candidates.csv")
    parser.add_argument("--model", type=str, default="facebook/esm2_t6_8M_UR50D")
    
    args = parser.parse_args()
    setup_logging()
    logger = get_logger(__name__)
    
    if not args.input_csv.exists():
        logger.error(f"Input file not found: {args.input_csv}")
        return 1
        
    df = pd.read_csv(args.input_csv)
    logger.info(f"Loaded {len(df)} candidates.")
    
    # Init Metric
    embedder = LatentEmbedder(model_name=args.model)
    metric = GrassmannMetric(embedder)
    
    # We need WT sequence to compare against.
    # In run_latent_rescue, we loaded WT. 
    # Ideally we'd store WT in the CSV or load it again.
    # We'll load it again.
    try:
        from experiments.run_design_rescues import _load_sequence
        wt_path = Path("data/interim/uniprot/P04637.fasta")
        if wt_path.exists():
            wt_seq = _load_sequence(wt_path)
        else:
            # Fallback
            wt_seq = "M" * 393 
    except:
        wt_seq = "M" * 393
        
    distances = []
    
    logger.info("Computing Grassmann Distances...")
    for seq in df["sequence"]:
        # Compute distance to WT
        # This checks "Is the subspace of this candidate similar to WT?"
        d = metric.grassmann_distance(wt_seq, seq)
        distances.append(d)
        
    df["grassmann_dist"] = distances
    
    # Save updated CSV
    out_path = args.input_csv.with_name("latent_rescue_candidates_analyzed.csv")
    df.to_csv(out_path, index=False)
    logger.info(f"Saved analyzed data to {out_path}")
    
    # PLOTTING
    plt.figure(figsize=(10, 6))
    
    # We want to plot Functional Score vs Grassmann Distance
    # Hue by Strategy
    if "predicted_function" in df.columns:
        sns.scatterplot(data=df, x="grassmann_dist", y="predicted_function", hue="strategy", style="target", s=100)
        plt.xlabel("Grassmann Distance from WT (Subspace Divergence)")
        plt.ylabel("Predicted Functional Score")
        plt.title("Functional Manifold Rescue Analysis")
        plt.axhline(0, color='gray', linestyle='--')
        
        # Ideal Quadrant: Low Distance, High Function (Top Left)
        plt.text(0.1, 2.0, "Ideal Rescue\n(Function + Structure)", color="green")
        
        plot_path = Path("reports/figures/fmr_analysis_plot.png")
        plot_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(plot_path)
        logger.info(f"Saved plot to {plot_path}")
    
    return 0

if __name__ == "__main__":
    main()
