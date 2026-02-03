"""
Analyze whether DNA contacts appear in Pareto front when constraint is disabled.
"""

import pandas as pd
import json
from pathlib import Path

# DNA-contacting positions (from 1TSR analysis)
DNA_CONTACTS = {
    # Strong contacts
    165, 167, 168, 249,
    # Medium contacts
    248, 250, 241, 120, 280, 282, 283, 273,
    # Weak contacts
    159, 277, 281, 286
}

def parse_mutations(rescue_str):
    """Extract positions from rescue mutations."""
    positions = set()
    for mut in rescue_str.split(','):
        # Extract position (e.g., "M133L" -> 133)
        pos_str = ''.join(c for c in mut if c.isdigit())
        if pos_str:
            positions.add(int(pos_str))
    return positions

def main():
    # Load candidates and Pareto front
    target = "R175H"
    base_path = Path("Data/processed/rescues") / target

    candidates_df = pd.read_parquet(base_path / "candidates.parquet")
    pareto_df = pd.read_parquet(base_path / "pareto.parquet")

    print(f"=== R175H Analysis (NO 8Å CONSTRAINT) ===\n")
    print(f"Total candidates: {len(candidates_df)}")
    print(f"Pareto front size: {len(pareto_df)}\n")

    # Check all candidates
    candidates_with_dna = []
    for idx, row in candidates_df.iterrows():
        positions = parse_mutations(row['rescue_mutations'])
        dna_positions = positions & DNA_CONTACTS
        if dna_positions:
            candidates_with_dna.append({
                'rescue': row['rescue_mutations'],
                'dna_positions': dna_positions,
                'risk': row['risk'],
                'ddg_gain': row['ddg_gain'],
                'in_pareto': row['rescue_mutations'] in pareto_df['rescue_mutations'].values
            })

    print(f"### Candidates with DNA contacts: {len(candidates_with_dna)} / {len(candidates_df)}")

    if candidates_with_dna:
        print("\nFirst 10 examples:")
        for i, item in enumerate(candidates_with_dna[:10], 1):
            pareto_flag = "⭐ IN PARETO" if item['in_pareto'] else ""
            print(f"  {i}. {item['rescue']:<20} DNA pos: {sorted(item['dna_positions'])} "
                  f"risk={item['risk']:.3f} ΔΔG={item['ddg_gain']:.2f} {pareto_flag}")

    # Check Pareto front specifically
    pareto_with_dna = []
    for idx, row in pareto_df.iterrows():
        positions = parse_mutations(row['rescue_mutations'])
        dna_positions = positions & DNA_CONTACTS
        if dna_positions:
            pareto_with_dna.append({
                'rescue': row['rescue_mutations'],
                'dna_positions': dna_positions,
                'risk': row['risk'],
                'ddg_gain': row['ddg_gain'],
                'rank': idx + 1  # Use index as rank since no rank column
            })

    print(f"\n### Pareto front with DNA contacts: {len(pareto_with_dna)} / {len(pareto_df)}")

    if pareto_with_dna:
        print("\n⚠️  DNA-CONTACTING RESCUES IN PARETO FRONT:")
        for item in pareto_with_dna:
            print(f"  Rank {item['rank']:3d}: {item['rescue']:<20} DNA pos: {sorted(item['dna_positions'])} "
                  f"risk={item['risk']:.3f} ΔΔG={item['ddg_gain']:.2f}")
    else:
        print("\n✅ NO DNA-contacting rescues in Pareto front")
        print("   → Algorithm naturally avoids DNA contacts (surface exposure + low stability)")

    # Check top 10
    print(f"\n### Top 10 Pareto rescues:")
    for idx, row in pareto_df.head(10).iterrows():
        positions = parse_mutations(row['rescue_mutations'])
        dna_flag = "⚠️  DNA CONTACT" if positions & DNA_CONTACTS else "✅"
        rank = idx + 1
        print(f"  {rank:3d}. {row['rescue_mutations']:<20} risk={row['risk']:.3f} ΔΔG={row['ddg_gain']:6.2f} {dna_flag}")

    # Summary
    print(f"\n=== CONCLUSION ===")
    if len(pareto_with_dna) == 0:
        print("✅ Pareto optimization NATURALLY avoids DNA contacts")
        print("   → The 8Å constraint is REDUNDANT (algorithm discovers this on its own)")
        print("   → Reason: DNA contacts are surface-exposed → low stability contribution")
    else:
        print(f"⚠️  Pareto front contains {len(pareto_with_dna)} DNA-contacting rescues")
        print("   → The 8Å constraint is NECESSARY (algorithm does not naturally avoid)")
        print("   → We were WRONG to claim emergent behavior")

    # Save results
    results = {
        'total_candidates': len(candidates_df),
        'pareto_size': len(pareto_df),
        'candidates_with_dna': len(candidates_with_dna),
        'pareto_with_dna': len(pareto_with_dna),
        'pareto_with_dna_list': [
            {
                'rescue': item['rescue'],
                'dna_positions': sorted(list(item['dna_positions'])),
                'rank': item['rank'],
                'risk': item['risk'],
                'ddg_gain': item['ddg_gain']
            }
            for item in pareto_with_dna
        ]
    }

    with open(base_path / "dna_contact_test_results.json", 'w') as f:
        json.dump(results, f, indent=2)

    print(f"\nResults saved to: {base_path / 'dna_contact_test_results.json'}")

if __name__ == "__main__":
    main()
