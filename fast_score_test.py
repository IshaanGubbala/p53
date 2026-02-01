#!/usr/bin/env python3

import sys
import time
from pathlib import Path

sys.path.insert(0, '/Users/ishaangubbala/Documents/p53/src')

from experiments.run_score_variants import _paths_from_config
from scoring.fast_runner import score_mutation_set
from scoring.common import canonicalize_mutation_set
import pandas as pd
from core.logging import get_logger
from core.config import load_configs

def run_fast_score():
    """Quick test of fast scoring"""
    print("Starting fast scoring test...")
    
    # Load configs
    configs = load_configs()
    paths = _paths_from_config(configs.get("paths", {}))
    
    # Load variants
    variant_path = paths["interim"] / "tp53_missense_normalized.parquet"
    df = pd.read_parquet(variant_path)
    mutation_set = canonicalize_mutation_set(df)
    
    # Test with first few mutations
    test_mutations = list(mutation_set)[:10]
    test_set = {mut: {'mutations': [(100, 'A', 'V')]} for mut in test_mutations}
    
    # Get PDB
    base_pdb = Path('/Volumes/Lexar/Data/raw/alphafold/AF-P04637-F1-model_v6.pdb')
    
    print(f"Scoring {len(test_mutations)} mutations...")
    
    # Score them
    start = time.time()
    results = score_mutation_set(test_set, base_pdb, Path('.'), {}, Path('.'))
    elapsed = time.time() - start
    
    print(f"Completed {len(test_mutations)} mutations in {elapsed:.2f}s")
    print(f"Average: {elapsed/len(test_mutations):.3f}s per mutation")
    
    return results

if __name__ == "__main__":
    run_fast_score()