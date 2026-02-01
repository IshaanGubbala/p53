"""
Test EvoEF2-based functional scoring on R175H top rescues.

This demonstrates the improvement from heuristic to energy-based scoring:
- Uses actual EvoEF2 ComputeBinding for DNA affinity and interface stability
- Same energy function across all dimensions (folding, binding, interface)
- Quantitative ΔΔG values in kcal/mol
"""

import pandas as pd
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.scoring.functional.functional_score_evoef2 import score_rescue_candidates_evoef2
from src.scoring.functional.evoef2_binding import load_evoef2_config


def main():
    print("=" * 80)
    print("EvoEF2-Based Functional Scoring Test: R175H Top Rescues")
    print("=" * 80)
    print()

    # Load R175H Pareto front
    pareto_path = Path("Data/processed/rescues/R175H/pareto.parquet")

    if not pareto_path.exists():
        print(f"Error: {pareto_path} not found")
        print("Please run rescue design first:")
        print("  python -m experiments.run_design_rescues --targets R175H")
        sys.exit(1)

    pareto_df = pd.read_parquet(pareto_path)

    print(f"Loaded {len(pareto_df)} Pareto-optimal rescues from stability-only optimization")
    print()
    print("=" * 80)
    print("Original Top 5 (ranked by ΔΔG_folding only)")
    print("=" * 80)

    # Show original top 5
    top5 = pareto_df.head(5)
    print(f"{'Rank':<6}{'Rescue':<20}{'ΔΔG_folding':<15}{'Risk':<10}")
    print("-" * 50)
    for idx, row in top5.iterrows():
        print(f"{idx+1:<6}{row['rescue_mutations']:<20}{row['ddg_gain']:<15.2f}{row['risk']:<10.3f}")

    print()
    print("=" * 80)
    print("Calculating EvoEF2-Based Functional Scores")
    print("=" * 80)
    print()
    print("For each rescue, computing:")
    print("  1. ΔΔG_folding  (monomer stability, already calculated)")
    print("  2. ΔΔG_binding  (DNA binding affinity, EvoEF2 ComputeBinding)")
    print("  3. ΔΔG_interface (tetramer stability, EvoEF2 ComputeBinding)")
    print("  4. Risk score   (MSA conservation, already calculated)")
    print()
    print("This will take 2-3 minutes per rescue (EvoEF2 builds mutants + calculates energies)")
    print("Using cached results when available...")
    print()

    # Load EvoEF2 config
    evoef2_cfg = load_evoef2_config()

    # Score top 5 only (full scoring takes ~15 min for all 162)
    scored_df = score_rescue_candidates_evoef2(
        top5,
        evoef2_cfg=evoef2_cfg,
        max_candidates=5,  # Limit to top 5 for testing
        verbose=True,
        cache_dir="Data/processed/cache/functional_scoring"
    )

    print()
    print("=" * 80)
    print("Results: EvoEF2-Based Functional Scores")
    print("=" * 80)
    print()

    # Sort by functional score
    scored_df = scored_df.sort_values('functional_score', ascending=False, na_position='last')

    # Show results table
    print(f"{'Rank':<6}{'Rescue':<15}{'Func':<8}{'Category':<12}{'ΔΔG_fold':<11}{'ΔΔG_bind':<11}{'ΔΔG_int':<11}{'Risk':<8}")
    print("-" * 90)

    for idx, (_, row) in enumerate(scored_df.iterrows(), 1):
        func_score = row.get('functional_score', float('nan'))
        cat = row.get('overall_category', 'N/A')
        ddg_fold = row.get('ddg_folding', row['ddg_gain'])
        ddg_bind = row.get('ddg_binding', float('nan'))
        ddg_int = row.get('ddg_interface', float('nan'))

        print(f"{idx:<6}"
              f"{row['rescue_mutations']:<15}"
              f"{func_score:<8.3f}"
              f"{cat:<12}"
              f"{ddg_fold:<11.2f}"
              f"{ddg_bind:<11.2f}"
              f"{ddg_int:<11.2f}"
              f"{row['risk']:<8.3f}")

    print()
    print("=" * 80)
    print("Analysis")
    print("=" * 80)
    print()

    # Count by category
    excellent = scored_df[scored_df['overall_category'] == 'excellent']
    good = scored_df[scored_df['overall_category'] == 'good']
    acceptable = scored_df[scored_df['overall_category'] == 'acceptable']
    poor = scored_df[scored_df['overall_category'] == 'poor']

    print(f"Category Distribution (of top 5 stability-optimized rescues):")
    print(f"  Excellent (all dimensions good):  {len(excellent)}")
    print(f"  Good (≥2 dimensions good):        {len(good)}")
    print(f"  Acceptable (≥1 dimension good):   {len(acceptable)}")
    print(f"  Poor (any dimension bad):         {len(poor)}")
    print()

    # Check for failures
    bad_binding = scored_df[scored_df['binding_category'] == 'bad']
    bad_interface = scored_df[scored_df['interface_category'] == 'bad']

    if len(bad_binding) > 0:
        print(f"⚠️  Rescues that FAIL DNA binding filter ({len(bad_binding)}):")
        for _, row in bad_binding.iterrows():
            print(f"  - {row['rescue_mutations']}: ΔΔG_binding = {row['ddg_binding']:.2f} kcal/mol")
        print()

    if len(bad_interface) > 0:
        print(f"⚠️  Rescues that FAIL tetramer interface filter ({len(bad_interface)}):")
        for _, row in bad_interface.iterrows():
            print(f"  - {row['rescue_mutations']}: ΔΔG_interface = {row['ddg_interface']:.2f} kcal/mol")
        print()

    if len(excellent) > 0:
        print(f"✅ Rescues rated EXCELLENT (all 3 dimensions good):")
        for _, row in excellent.iterrows():
            print(f"  ⭐ {row['rescue_mutations']}: "
                  f"functional_score = {row['functional_score']:.3f}")
        print()

    # Ranking changes
    print("=" * 80)
    print("Ranking Changes: Stability-Only vs Functional")
    print("=" * 80)
    print()

    for idx, (_, row) in enumerate(scored_df.iterrows(), 1):
        rescue = row['rescue_mutations']
        old_rank = top5[top5['rescue_mutations'] == rescue].index[0] + 1
        new_rank = idx
        change = old_rank - new_rank

        arrow = "↑" if change > 0 else "↓" if change < 0 else "→"
        symbol = "✅" if row.get('overall_category') in ['excellent', 'good'] else "⚠️"

        print(f"  {symbol} {rescue:<20} Old: {old_rank}  {arrow}  New: {new_rank}  "
              f"({change:+2d})  [{row.get('overall_category', 'N/A')}]")

    # Save results
    output_path = Path("Data/processed/functional_scores") / "R175H_top5_evoef2.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    scored_df.to_csv(output_path, index=False)

    print()
    print("=" * 80)
    print(f"✅ Results saved to: {output_path}")
    print("=" * 80)
    print()

    # Key insights
    print("Key Insights:")
    print()
    print("1. Stability-only optimization may promote rescues that:")
    print("   - Stabilize the monomer BUT disrupt DNA binding")
    print("   - Stabilize the monomer BUT destabilize the tetramer")
    print()
    print("2. Functional scoring filters these false positives by:")
    print("   - Using EvoEF2 ComputeBinding to calculate ΔΔG_binding")
    print("   - Using EvoEF2 ComputeBinding to calculate ΔΔG_interface")
    print("   - Ensuring consistency (same energy function for all dimensions)")
    print()
    print("3. 'Excellent' rescues pass all 3 dimensions:")
    print("   - ΔΔG_folding ≤ 0 (stabilizing)")
    print("   - ΔΔG_binding ≤ 0 (preserves/improves DNA binding)")
    print("   - ΔΔG_interface ≤ 0 (preserves/improves tetramer)")
    print()


if __name__ == "__main__":
    main()
