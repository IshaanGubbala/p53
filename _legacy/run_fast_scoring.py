#!/usr/bin/env python3

import pandas as pd
import time
from pathlib import Path
import sys
import os

# Add to path for imports
sys.path.insert(0, '/Users/ishaangubbala/Documents/p53/src')

from scoring.fast_runner import score_mutation_set
from scoring.common import canonicalize_mutation_set

def main():
    print('=== FAST SCORING FOR TP53 MUTATIONS ===')
    
    # Load all variants
    try:
        df = pd.read_parquet('/Volumes/Lexar/Data/interim/tp53_missense_normalized.parquet')
    except Exception as e:
        print(f'Error loading variants: {e}')
        return
    
    print(f'Loaded {len(df)} variants from parquet')
    
    # Create mutation list from dataframe
    mutations = []
    for _, row in df.iterrows():
        if pd.notna(row['pos']) and pd.notna(row['ref']) and pd.notna(row['alt']):
            mut_str = f"{row['ref']}{int(row['pos'])}{row['alt']}"
            mutations.append(mut_str)
    
    print(f'Converted {len(mutations)} mutations to string format')
    
    # Score in batches
    batch_size = 100
    base_pdb = Path('/Volumes/Lexar/Data/raw/alphafold/AF-P04637-F1-model_v6.pdb')
    
    if not base_pdb.exists():
        print(f'PDB file not found: {base_pdb}')
        return
    
    all_results = {}
    total_start = time.time()
    
    for i in range(0, len(mutations), batch_size):
        batch_mutations = mutations[i:i + batch_size]
        
        # Create simple test set
        test_set = {}
        for mut in batch_mutations:
            pos = int(''.join(filter(str.isdigit, mut)))
            ref = ''.join(filter(str.isalpha, mut))[0]
            alt = ''.join(filter(str.isalpha, mut))[-1]
            test_set[mut] = {'mutations': [(pos, ref, alt)]}
        
        try:
            results = score_mutation_set(test_set, base_pdb, Path('.'), {}, Path('.'))
            all_results.update(results)
            
            batch_elapsed = time.time() - total_start
            progress = (i + len(batch_mutations)) / len(mutations) * 100
            print(f'Progress: {progress:.1f}% ({i + len(batch_mutations)}/{len(mutations)}) - Batch {i//batch_size + 1}: {len(batch_mutations)} mutations in {batch_elapsed:.2f}s')
            
        except Exception as e:
            print(f'Error scoring batch {i//batch_size + 1}: {e}')
    
    total_elapsed = time.time() - total_start
    print(f'\\n=== COMPLETED {len(all_results)} mutations ===')
    print(f'Total time: {total_elapsed:.2f}s')
    print(f'Average speed: {total_elapsed/len(all_results):.3f}s per mutation')
    
    # Save results
    if all_results:
        import json
        output_path = Path('/Volumes/Lexar/Data/processed/variant_scores_fast.json')
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        with output_path.open('w') as f:
            json.dump(all_results, f, indent=2)
        
        print(f'Results saved to: {output_path}')

if __name__ == "__main__":
    main()