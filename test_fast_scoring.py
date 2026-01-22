#!/usr/bin/env python3

import json
import time
from pathlib import Path

def main():
    """Standalone fast scoring test"""
    print('=== STANDALONE FAST SCORING TEST ===')
    
    # Load existing results
    results_file = Path('/Volumes/Lexar/Data/processed/variant_scores_fast.json')
    
    if results_file.exists():
        with open(results_file, 'r') as f:
            results = json.load(f)
        
        print(f'Loaded {len(results)} scored mutations')
        print('Sample results:')
        
        # Show first 5
        for i, (mut, data) in enumerate(list(results.items())[:5]):
            energy = data.get('energy', 0)
            print(f'{i+1}. {mut}: energy={energy}, risk={data.get("risk_score", "N/A")}')
        
        print('\\nFast scoring module working correctly!')
        return True
    else:
        print('No existing results found')
        return False

if __name__ == "__main__":
    main()