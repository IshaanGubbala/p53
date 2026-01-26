"""
Run DNA binding affinity analysis on rescue candidates.

This analysis identifies rescues that might destroy DNA binding
despite stabilizing the protein - a critical validation step.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import yaml

from src.core.logging import get_logger
from src.eval.dna_binding_analysis import (
    analyze_rescue_dna_binding,
    identify_dna_contacts,
    classify_dna_contact_strength,
    summarize_dna_binding_analysis,
)

logger = get_logger(__name__)


def run(args, configs: dict) -> int:
    """Run DNA binding analysis on rescue candidates."""
    logger.info("Starting DNA binding affinity analysis...")

    # Paths
    paths_cfg = configs.get("paths", {})
    data_cfg = paths_cfg.get("data", {})
    processed_dir = Path(data_cfg.get("processed", "Data/processed"))
    raw_dir = Path(data_cfg.get("raw", "Data/raw"))
    reports_dir = Path(paths_cfg.get("reports_dir", "reports"))

    rescues_root = processed_dir / "rescues"
    dna_complex_pdb = raw_dir / "experimental_pdbs" / "1TSR.pdb"

    if not dna_complex_pdb.exists():
        logger.error(f"DNA-bound structure not found: {dna_complex_pdb}")
        logger.info("Please download 1TSR.pdb (p53-DNA complex) to Data/raw/experimental_pdbs/")
        return 1

    if not rescues_root.exists():
        logger.error(f"Rescues directory not found: {rescues_root}")
        return 1

    # Output directory
    dna_binding_dir = reports_dir / "dna_binding"
    dna_binding_dir.mkdir(parents=True, exist_ok=True)

    # First, identify DNA contacts in the structure
    logger.info(f"Analyzing DNA contacts in {dna_complex_pdb.name}...")
    dna_contacts = identify_dna_contacts(
        dna_complex_pdb,
        protein_chains=["A"],  # Use chain A
        dna_chains=["E", "F"],  # DNA chains
    )

    # Classify contact strengths
    contact_classifications = classify_dna_contact_strength(dna_contacts)

    # Save DNA contact information
    import json
    contacts_json = dna_binding_dir / "dna_contacts.json"
    with open(contacts_json, "w") as f:
        json.dump({
            "contacts": dna_contacts,
            "classifications": contact_classifications,
        }, f, indent=2)
    logger.info(f"Saved DNA contact analysis to {contacts_json}")

    # Print DNA contact summary
    print(f"\n{'='*80}")
    print(f"DNA-Contacting Residues (from 1TSR structure)")
    print(f"{'='*80}")
    print(f"{'Position':<10} {'Residue':<10} {'Strength':<12} {'Contacts':<10} {'Min Dist (Å)'}")
    print(f"{'-'*80}")

    for pos in sorted(dna_contacts.keys()):
        res_name = dna_contacts[pos]["residue_name"]
        strength = contact_classifications.get(pos, "unknown")
        n_contacts = dna_contacts[pos]["n_contacts"]
        min_dist = dna_contacts[pos]["min_distance"]
        print(f"{pos:<10} {res_name:<10} {strength:<12} {n_contacts:<10} {min_dist:.2f}")

    # Process each target
    all_results = []

    for target_dir in sorted(rescues_root.iterdir()):
        if not target_dir.is_dir():
            continue

        target = target_dir.name
        pareto_path = target_dir / "pareto.parquet"

        if not pareto_path.exists():
            logger.warning(f"No pareto file for {target}, skipping")
            continue

        logger.info(f"\nAnalyzing DNA binding impact for {target}...")

        # Load rescue candidates
        pareto_df = pd.read_parquet(pareto_path)

        # Sort by Pareto status and quality
        pareto_df = pareto_df.sort_values(
            ["is_pareto", "n_rescue", "ddg_gain", "risk"],
            ascending=[False, True, True, True]
        )

        try:
            # Analyze DNA binding impact
            dna_binding_df = analyze_rescue_dna_binding(
                pareto_df.head(50),  # Analyze top 50
                dna_complex_pdb=dna_complex_pdb,
                protein_chains=["A"],
                dna_chains=["E", "F"],
                filter_high_risk=False,  # Don't filter yet, just annotate
            )

            # Save results
            output_csv = dna_binding_dir / f"{target}_dna_binding.csv"
            dna_binding_df.to_csv(output_csv, index=False)
            logger.info(f"Saved DNA binding results to {output_csv}")

            # Generate summary
            summary = summarize_dna_binding_analysis(
                dna_binding_df,
                output_path=dna_binding_dir / f"{target}_dna_binding_summary.json"
            )

            # Print summary
            print(f"\n{'='*80}")
            print(f"DNA Binding Analysis: {target}")
            print(f"{'='*80}")
            print(f"  Rescues analyzed: {summary['n_rescues_analyzed']}")
            print(f"  Rescues contacting DNA: {summary['n_contact_dna']}")
            print(f"  Mean impact score: {summary['mean_impact_score']:.3f}")
            print(f"  Max impact score: {summary['max_impact_score']:.3f}")

            print(f"\n  Risk Distribution:")
            for risk_level, count in sorted(summary['risk_distribution'].items()):
                print(f"    {risk_level:12s}: {count}")

            # Identify dangerous rescues (high DNA binding impact)
            dangerous = dna_binding_df[
                dna_binding_df["dna_binding_risk"].isin(["critical", "high"])
            ]

            if not dangerous.empty:
                print(f"\n  ⚠️  WARNING: {len(dangerous)} rescues with high DNA binding risk:")
                for idx, row in dangerous.head(5).iterrows():
                    print(f"    - {row['rescue_mutations']:20s} | "
                          f"Impact: {row['dna_binding_impact_score']:.2f} | "
                          f"Risk: {row['dna_binding_risk']}")

            # Identify safe rescues
            safe = dna_binding_df[
                dna_binding_df["dna_binding_risk"].isin(["safe", "low"])
            ]

            if not safe.empty:
                print(f"\n  ✅ Top 3 safe rescues (low DNA binding impact):")
                for idx, row in safe.head(3).iterrows():
                    print(f"    - {row['rescue_mutations']:20s} | "
                          f"Impact: {row['dna_binding_impact_score']:.2f} | "
                          f"Risk: {row['dna_binding_risk']}")

            all_results.append({
                "target": target,
                "summary": summary,
                "n_dangerous": len(dangerous),
                "n_safe": len(safe),
            })

        except Exception as e:
            logger.error(f"Failed to analyze {target}: {e}", exc_info=True)
            continue

    # Generate consolidated report
    if all_results:
        print(f"\n{'='*80}")
        print("Consolidated DNA Binding Summary")
        print(f"{'='*80}")

        for result in all_results:
            target = result["target"]
            summary = result["summary"]
            print(f"\n{target}:")
            print(f"  DNA-contacting rescues: {summary['n_contact_dna']}")
            print(f"  High-risk rescues: {result['n_dangerous']}")
            print(f"  Safe rescues: {result['n_safe']}")
            print(f"  Mean impact score: {summary['mean_impact_score']:.3f}")

        # Save consolidated CSV
        all_dna_binding = []
        for target_dir in sorted(rescues_root.iterdir()):
            if not target_dir.is_dir():
                continue
            target = target_dir.name
            csv_path = dna_binding_dir / f"{target}_dna_binding.csv"
            if csv_path.exists():
                df = pd.read_csv(csv_path)
                df["target"] = target
                all_dna_binding.append(df)

        if all_dna_binding:
            consolidated_df = pd.concat(all_dna_binding, ignore_index=True)
            consolidated_csv = dna_binding_dir / "all_targets_dna_binding.csv"
            consolidated_df.to_csv(consolidated_csv, index=False)
            logger.info(f"Saved consolidated results to {consolidated_csv}")

    print(f"\n{'='*80}")
    print(f"DNA binding analysis complete!")
    print(f"Results saved to: {dna_binding_dir}")
    print(f"{'='*80}\n")

    return 0


if __name__ == "__main__":
    # Load configs
    config_files = [
        "configs/p53.yaml",
        "configs/paths.yaml",
    ]

    configs = {}
    for cfg_file in config_files:
        cfg_path = Path(cfg_file)
        if cfg_path.exists():
            with open(cfg_path) as f:
                cfg_data = yaml.safe_load(f) or {}
                configs.update(cfg_data)

    sys.exit(run(None, configs))
