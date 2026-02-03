
import pandas as pd
from pathlib import Path
import yaml
import sys
import torch

# Add project root to path
sys.path.append(str(Path.cwd()))

from src.scoring.rasp_scorer import score_mutation_set_rasp
from src.core.logging import get_logger

def main():
    logger = get_logger(__name__)
    summary_path = Path("reports/executive_summary.csv")
    if not summary_path.exists():
        logger.error(f"Summary not found: {summary_path}")
        return

    df = pd.read_csv(summary_path)
    
    # Load scoring config to get rasp settings and base PDB
    with open("configs/scoring.yaml", "r") as f:
        scoring_cfg = yaml.safe_load(f)
    
    rasp_cfg = scoring_cfg.get("rasp", {})
    if not rasp_cfg.get("enabled", False):
        logger.warning("RaSP is disabled in config. Enabling for this run.")
        rasp_cfg["enabled"] = True
    
    # Base PDB from config (repaired AF model)
    base_pdb = Path(scoring_cfg["evoef2"]["repaired_pdb"])
    cache_dir = Path("Data/processed/cache")
    
    logger.info(f"Running RaSP validation on {len(df)} candidates using {base_pdb}")
    
    rasp_scores = []
    for idx, row in df.iterrows():
        mutations = [row['target']] # The seed mutation
        if pd.notna(row['rescue_mutations']):
            rescues = [m.strip() for m in str(row['rescue_mutations']).split(",")]
            mutations.extend(rescues)
        
        logger.info(f"Scoring set {idx+1}/{len(df)}: {mutations}")
        try:
            # Note: score_mutation_set_rasp returns total ddg for the set
            score = score_mutation_set_rasp(mutations, base_pdb, cache_dir, rasp_cfg)
            rasp_scores.append(round(score, 3))
        except Exception as e:
            logger.error(f"Failed to score {mutations}: {e}")
            rasp_scores.append(None)
            
    df['RaSP_ddg'] = rasp_scores
    
    # Save back
    df.to_csv(summary_path, index=False)
    logger.info(f"Updated {summary_path} with RaSP scores")
    print(df[['target', 'rescue_mutations', 'ddg_seed', 'ddg_gain', 'RaSP_ddg']])

if __name__ == "__main__":
    main()
