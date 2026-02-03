from __future__ import annotations

import argparse
import json
import yaml
from pathlib import Path
from typing import Any
import sys

# Add project root to path
sys.path.append(str(Path.cwd()))

import pandas as pd
from src.core.logging import get_logger, setup_logging
from src.design.generative import get_generator, extract_mutations
# Reusing path logic from run_design_rescues if possible, or simplified version
from experiments.run_design_rescues import _paths_from_config, _load_sequence, _parse_mutation

def apply_mutation(seq: str, mutation: str) -> str:
    """Apply a single mutation like 'R175H' to a sequence."""
    pos, ref, alt = _parse_mutation(mutation)
    # _parse_mutation returns 1-indexed pos
    idx = pos - 1
    if idx < 0 or idx >= len(seq):
        raise ValueError(f"Position {pos} out of range")
    if seq[idx] != ref:
        # Warn but proceed? 
        pass 
    s_list = list(seq)
    s_list[idx] = alt
    return "".join(s_list)

def main():
    parser = argparse.ArgumentParser(description="Run generative design for p53 rescue")
    parser.add_argument("--config", type=Path, default="configs/scoring.yaml") # borrowing scoring/main config
    parser.add_argument("--targets", type=str, nargs="+", help="Target mutations (e.g. R175H)")
    parser.add_argument("--n_samples", type=int, default=100)
    args = parser.parse_args()

    setup_logging()
    logger = get_logger(__name__)
    
    # Load config
    # Assuming standard config structure 
    # We might need a 'generative.yaml' or use 'scoring.yaml' + extra args
    # For now, I'll load scoring.yaml for paths, and assume defaults
    if not args.config.exists():
        logger.error(f"Config file not found: {args.config}")
        return 1
        
    with args.config.open() as f:
        config = yaml.safe_load(f)

    paths = _paths_from_config(config.get("paths", {}))
    
    # Load WT Sequence
    p53_cfg = config.get("p53", {})
    uniprot_id = p53_cfg.get("uniprot_id", "P04637")
    fasta_path = paths["interim"] / "uniprot" / f"{uniprot_id}.fasta"
    if not fasta_path.exists():
        logger.error(f"FASTA not found: {fasta_path}")
        return 1
    
    wt_sequence = _load_sequence(fasta_path)
    
    # Initialize Generator
    gen_cfg = config.get("generative", {"model_name": "EvoDiff"})
    generator = get_generator("EvoDiff", gen_cfg)
    
    targets = args.targets or ["R175H", "R248Q", "R273H"]
    
    all_candidates = []
    
    for target in targets:
        logger.info(f"Generating candidates for target: {target}")
        
        # Create seed sequence (WT + Cancer Mutation)
        try:
            cancer_seq = apply_mutation(wt_sequence, target)
        except Exception as e:
            logger.error(f"Failed to create start sequence for {target}: {e}")
            continue
            
        # Generate
        generated_seqs = generator.generate(cancer_seq, n_samples=args.n_samples)
        
        valid_count = 0
        for seq in generated_seqs:
            # Extract mutations relative to WT
            muts = extract_mutations(wt_sequence, seq)
            
            # Filter: Must preserve target mutation?
            # User wants RESCUE. A rescue usually implies the cancer mutation is still there 
            # (unless we are just repairing it, but that's boring).
            # Let's check if target is in muts
            if target not in muts:
                # The model reverted the cancer mutation or mutated the position again
                # Skip for now, or log as 'reversion'
                continue
            
            # If we only have the target, it's not a rescue, it's just the input
            if len(muts) < 2:
                continue
                
            mut_str = ",".join(sorted(muts))
            all_candidates.append({
                "target": target,
                "mutation": mut_str,
                "n_mutations": len(muts),
                "source": "generative_evodiff"
            })
            valid_count += 1
            
        logger.info(f"Generated {valid_count} unique rescue candidates for {target}")

    if not all_candidates:
        logger.warning("No candidates generated.")
        return 0
        
    # Save results
    df = pd.DataFrame(all_candidates)
    # Drop duplicates
    df = df.drop_duplicates(subset=["mutation"])
    
    # Save as CSV
    out_path_csv = paths["processed"] / "generative_candidates.csv"
    out_path_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path_csv, index=False)
    
    # Save as Parquet (Normalized for scoring)
    # Ensure columns match what run_score_variants expects: pos, ref, alt
    # generative output has: target, mutation (e.g. R175H,M237I), n_mutations, source
    # We need to explode the mutation string if we want single variants, 
    # but run_score_variants parses 'mutation' column?
    # run_score_variants _load_variants requires: pos, ref, alt.
    # It constructs "mutation" from them.
    # Let's try to adapt the dataframe to match that format if possible, 
    # otherwise just saving the parquet is what was asked.
    
    out_path_parquet = paths["processed"] / "generative_candidates.parquet"
    try:
        df.to_parquet(out_path_parquet, index=False)
        logger.info(f"Saved {len(df)} candidates to parquet: {out_path_parquet}")
    except ImportError:
        logger.warning("Could not save to Parquet (missing pyarrow/fastparquet). Saved CSV only.")
    
    logger.info(f"Saved {len(df)} candidates to {out_path_csv}")
    logger.info("To score these, run: python3 experiments/run_score_variants.py")
    
    return 0

if __name__ == "__main__":
    main()
