"""
Test functional scoring on real R175H rescue candidates.

This script:
1. Loads the R175H Pareto front (top rescues)
2. Applies functional scoring (DNA binding + interface + folding + risk)
3. Compares old rankings (stability-only) vs new rankings (functional)
4. Identifies rescues that fail DNA/interface filters
"""

import pandas as pd
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.scoring.functional import score_rescue_candidates


def main():
    print("=== Functional Scoring Test: R175H Top Rescues ===\n")

    # Load R175H Pareto front
    pareto_path = Path("Data/processed/rescues/R175H/pareto.parquet")

    if not pareto_path.exists():
        print(f"Error: {pareto_path} not found")
        print("Please run rescue design first")
        sys.exit(1)

    pareto_df = pd.read_parquet(pareto_path)

    print(f"Loaded {len(pareto_df)} Pareto-optimal rescues")
    print(f"\n=== Original Top 10 (by stability) ===")

    # Show original top 10
    top10 = pareto_df.head(10)
    for idx, row in top10.iterrows():
        print(f"  {idx+1:2d}. {row['rescue_mutations']:<20} "
              f"ΔΔG={row['ddg_gain']:6.2f}  risk={row['risk']:.3f}")

    print(f"\n=== Calculating Functional Scores ===")
    print("This evaluates DNA binding + interface + folding + risk...")
    print("(May take a few minutes for {len(top10)} rescues)\n")

    # Score top 10 only (for speed)
    scored_df = score_rescue_candidates(
        top10,
        verbose=True,
        cache_dir="Data/processed/cache/functional_scoring"
    )

    print(f"\n=== Results ===\n")

    # Sort by functional score
    scored_df = scored_df.sort_values('functional_score', ascending=False)

    # Show results
    print(f"{'Rank':<6}{'Rescue':<20}{'Func Score':<12}{'Category':<12}"
          f"{'ΔΔG_fold':<10}{'ΔΔG_bind':<10}{'ΔΔG_int':<10}{'Risk':<8}")
    print("-" * 95)

    for idx, (_, row) in enumerate(scored_df.iterrows(), 1):
        print(f"{idx:<6}"
              f"{row['rescue_mutations']:<20}"
              f"{row.get('functional_score', np.nan):<12.3f}"
              f"{row.get('overall_category', 'N/A'):<12}"
              f"{row.get('ddg_folding', row['ddg_gain']):<10.2f}"
              f"{row.get('ddg_binding', np.nan):<10.2f}"
              f"{row.get('ddg_interface', np.nan):<10.2f}"
              f"{row['risk']:<8.3f}")

    # Identify failures
    print(f"\n=== Filtering Analysis ===\n")

    bad_binding = scored_df[scored_df['binding_category'] == 'bad']
    bad_interface = scored_df[scored_df['interface_category'] == 'bad']

    print(f"Rescues that FAIL DNA binding filter: {len(bad_binding)}")
    if len(bad_binding) > 0:
        for _, row in bad_binding.iterrows():
            print(f"  - {row['rescue_mutations']}: ΔΔG_binding = {row['ddg_binding']:.2f}")

    print(f"\nRescues that FAIL interface filter: {len(bad_interface)}")
    if len(bad_interface) > 0:
        for _, row in bad_interface.iterrows():
            print(f"  - {row['rescue_mutations']}: ΔΔG_interface = {row['ddg_interface']:.2f}")

    # Excellent rescues
    excellent = scored_df[scored_df['overall_category'] == 'excellent']
    print(f"\nRescues rated 'EXCELLENT' (all dimensions good): {len(excellent)}")
    if len(excellent) > 0:
        for _, row in excellent.iterrows():
            print(f"  ⭐ {row['rescue_mutations']}: "
                  f"Func={row['functional_score']:.3f}")

    # Compare rankings
    print(f"\n=== Ranking Changes ===\n")
    print("Comparing stability-only vs functional score rankings:\n")

    for idx, (_, row) in enumerate(scored_df.iterrows(), 1):
        rescue = row['rescue_mutations']
        old_rank = top10[top10['rescue_mutations'] == rescue].index[0] + 1
        new_rank = idx
        change = old_rank - new_rank

        arrow = "↑" if change > 0 else "↓" if change < 0 else "→"
        print(f"  {rescue:<20} {old_rank:2d} {arrow} {new_rank:2d}  "
              f"(change: {change:+3d})")

    # Save results
    output_path = Path("Data/processed/functional_scores") / "R175H_top10_functional.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    scored_df.to_csv(output_path, index=False)
    print(f"\n✅ Results saved to: {output_path}")


if __name__ == "__main__":
    import numpy as np
    main()
