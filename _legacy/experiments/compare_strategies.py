
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import json
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

import argparse

def load_search_results(processed_dir: Path, target: str):
    target_dir = processed_dir / "rescues" / target
    candidates_path = target_dir / "candidates.parquet"
    
    if not candidates_path.exists():
        logger.warning(f"Search results not found at {candidates_path}")
        return None
        
    df = pd.read_parquet(candidates_path)
    df["strategy"] = "search"
    df["target"] = target
    return df

def load_generative_results(processed_dir: Path, target: str):
    # check for scored version first
    gen_path = processed_dir / "generative_candidates_scored.parquet"
    if not gen_path.exists():
        gen_path = processed_dir / "generative_candidates.parquet"
    if not gen_path.exists():
        gen_path = processed_dir / "generative_candidates.csv"
    
    if not gen_path.exists():
        logger.warning(f"Generative results not found at {gen_path}")
        return None
    
    if gen_path.suffix == ".parquet":
        df = pd.read_parquet(gen_path)
    else:
        df = pd.read_csv(gen_path)
        
    # Filter by target if the column exists
    if "target" in df.columns:
        df = df[df["target"] == target].copy()
        
    df["strategy"] = "generative"
    return df

def plot_comparison(search_df, gen_df, target):
    combined = pd.concat([search_df, gen_df])
    
    # Filter for candidates that actually have scores
    # Search candidates have ddg and rasp_ddg
    # Gen candidates might have them if scored
    
    plt.figure(figsize=(10, 6))
    sns.scatterplot(data=combined, x="ddg", y="rasp_ddg", hue="strategy", style="strategy", alpha=0.6)
    
    plt.title(f"Search vs Generative Strategy - {target}")
    plt.xlabel("EvoEF2 ΔΔG (Stability)")
    plt.ylabel("RaSP ΔΔG (Deep Learning)")
    plt.axvline(x=0, color='gray', linestyle='--')
    plt.axhline(y=0, color='gray', linestyle='--')
    
    plt.tight_layout()
    out_path = f"reports/strategy_comparison_{target}.png"
    Path("reports").mkdir(exist_ok=True)
    plt.savefig(out_path)
    logger.info(f"Plot saved to {out_path}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", default="R175H")
    args = parser.parse_args()
    
    processed_dir = Path("Data/processed")
    
    search_df = load_search_results(processed_dir, args.target)
    gen_df = load_generative_results(processed_dir, args.target)
    
    if search_df is None and gen_df is None:
        logger.error("No results found for comparison.")
        return
        
    if search_df is not None:
        logger.info(f"Loaded {len(search_df)} search candidates for {args.target}")
    if gen_df is not None:
        logger.info(f"Loaded {len(gen_df)} generative candidates for {args.target}")

    if search_df is not None and gen_df is not None:
        # Check if gen_df has scores
        if "ddg" in gen_df.columns:
            plot_comparison(search_df, gen_df, args.target)
        else:
            logger.warning("Generative candidates not yet scored. Run: python3 experiments/run_score_variants.py --variants Data/processed/generative_candidates.parquet")


if __name__ == "__main__":
    main()
