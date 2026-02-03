#!/usr/bin/env python3
"""Create tiered rescue recommendations with mechanistic rationale.

Selects best single, double, and triple for each target based on
Pareto optimality and risk/gain tradeoffs.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

# Add repo root to Python path
repo_root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(repo_root))


def load_structure_data(interim_dir: Path) -> tuple[dict, dict]:
    """Load burial and conservation data for mechanistic rationale."""
    structure_dir = interim_dir / "structure"
    burial_path = structure_dir / "P04637_core_94_312_burial.json"

    burial_map = {}
    if burial_path.exists():
        burial_raw = json.loads(burial_path.read_text(encoding="utf-8"))
        # Convert categorical burial to numeric scores
        burial_scores = {"buried": 1.0, "partial": 0.5, "exposed": 0.0}
        burial_map = {int(k): burial_scores.get(v, 0.0) for k, v in burial_raw.items()}

    # Load conservation if available
    cons_path = interim_dir / "conservation.json"
    cons_map = {}
    if cons_path.exists():
        cons_map = json.loads(cons_path.read_text(encoding="utf-8"))
        cons_map = {int(k): float(v) for k, v in cons_map.items()}

    return burial_map, cons_map


def parse_mutation(mut_str: str) -> tuple[str, int, str]:
    """Parse mutation string like 'S95A' into (ref, pos, alt)."""
    import re
    match = re.match(r'([A-Z*])(\d+)([A-Z*])', mut_str)
    if not match:
        raise ValueError(f"Invalid mutation: {mut_str}")
    ref, pos, alt = match.groups()
    return ref, int(pos), alt


def get_mutation_rationale(mutation: str, burial_map: dict, cons_map: dict) -> str:
    """Generate mechanistic rationale for a mutation."""
    ref, pos, alt = parse_mutation(mutation)

    rationale_parts = []

    # Check burial
    burial = burial_map.get(pos, 0.0)
    if burial > 0.7:
        rationale_parts.append(f"buried core packing (burial={burial:.2f})")
    elif burial > 0.4:
        rationale_parts.append(f"semi-buried stabilization (burial={burial:.2f})")
    else:
        rationale_parts.append(f"surface stabilization (burial={burial:.2f})")

    # Check conservation
    if cons_map:
        cons = cons_map.get(pos, 0.0)
        if cons < 0.3:
            rationale_parts.append("low conservation (safe to mutate)")

    # Add mutation-specific rationale based on amino acid properties
    hydrophobic = {'A', 'V', 'L', 'I', 'M', 'F', 'W', 'P'}
    aromatic = {'F', 'Y', 'W', 'H'}
    small = {'A', 'G', 'S'}

    if ref in hydrophobic and alt == 'A':
        rationale_parts.append("alanine substitution for packing optimization")
    elif ref == 'S' and alt == 'A':
        rationale_parts.append("serine→alanine removes potential instability")
    elif ref == 'M' and alt == 'L':
        rationale_parts.append("Met→Leu: similar hydrophobic, more stable")
    elif ref in aromatic and alt in aromatic:
        rationale_parts.append("aromatic substitution maintains interactions")
    elif ref == 'R' and alt in {'K', 'Q', 'H'}:
        rationale_parts.append("conservative charged/polar substitution")
    elif ref == 'A' and alt in {'S', 'G'}:
        rationale_parts.append("small residue substitution for flexibility")

    return "; ".join(rationale_parts)


def select_best_by_complexity(df: pd.DataFrame, burial_map: dict, cons_map: dict) -> dict:
    """Select best single, double, and triple from Pareto front."""
    # First try Pareto front
    pareto = df[df["is_pareto"]].copy()

    # If not enough Pareto, fall back to all candidates
    if len(pareto) < 3:
        pareto = df.copy()

    results = {}

    for n in [1, 2, 3]:
        subset = pareto[pareto["n_rescue"] == n]

        if len(subset) == 0:
            continue

        # Select best: lowest risk, then best ddg_gain
        best = subset.sort_values(["risk", "ddg_gain"]).iloc[0]

        # Parse rescue mutations
        rescue_list = best["rescue_mutations"].split(",") if best["rescue_mutations"] else []

        # Generate rationale for each mutation
        rationales = {}
        for mut in rescue_list:
            rationales[mut] = get_mutation_rationale(mut, burial_map, cons_map)

        results[f"best_{n}mut"] = {
            "rescue_mutations": rescue_list,
            "full_mutations": best["full_mutations"],
            "ddg_seed": float(best["ddg_seed"]),
            "ddg_total": float(best["ddg_total"]),
            "ddg_gain": float(best["ddg_gain"]),
            "risk": float(best["risk"]),
            "risk_components": json.loads(best["risk_components"]),
            "is_pareto": bool(best["is_pareto"]),
            "rationale": rationales,
        }

    return results


def main() -> None:
    data_dir = repo_root / "Data"
    interim_dir = data_dir / "interim"
    processed_dir = data_dir / "processed"
    reports_dir = repo_root / "reports"
    tables_dir = reports_dir / "tables"

    # Load structure data for rationale
    burial_map, cons_map = load_structure_data(interim_dir)

    targets = ["R175H", "R248Q", "R273H"]

    print("=== Creating Tiered Rescue Recommendations ===\n")

    all_recommendations = {}

    for target in targets:
        print(f"Processing {target}...")

        target_dir = processed_dir / "rescues" / target
        candidates_path = target_dir / "candidates.parquet"

        if not candidates_path.exists():
            print(f"  Skipping {target}: no candidates found")
            continue

        df = pd.read_parquet(candidates_path)
        recommendations = select_best_by_complexity(df, burial_map, cons_map)
        all_recommendations[target] = recommendations

        # Print summary
        for key, rec in recommendations.items():
            n_muts = len(rec["rescue_mutations"])
            print(f"  {key}: {', '.join(rec['rescue_mutations'])}")
            print(f"    ΔΔG gain: {rec['ddg_gain']:.2f}, Risk: {rec['risk']:.3f}, Pareto: {rec['is_pareto']}")

        # Save individual target file
        target_output = tables_dir / f"top3_by_complexity_{target}.csv"

        rows = []
        for key, rec in recommendations.items():
            rows.append({
                "tier": key,
                "rescue_mutations": ",".join(rec["rescue_mutations"]),
                "full_mutations": rec["full_mutations"],
                "n_rescue": len(rec["rescue_mutations"]),
                "ddg_seed": rec["ddg_seed"],
                "ddg_total": rec["ddg_total"],
                "ddg_gain": rec["ddg_gain"],
                "risk": rec["risk"],
                "is_pareto": rec["is_pareto"],
                "rationale": json.dumps(rec["rationale"]),
            })

        pd.DataFrame(rows).to_csv(target_output, index=False)
        print(f"  Saved to {target_output}\n")

    # Save comprehensive JSON
    json_output = tables_dir / "tiered_recommendations.json"
    json_output.write_text(json.dumps(all_recommendations, indent=2), encoding="utf-8")
    print(f"✓ All recommendations saved to {json_output}")


if __name__ == "__main__":
    main()
