"""
Run tetramer interface analysis on rescue candidates.

This analysis identifies rescues that might disrupt p53 tetramerization,
causing dominant-negative effects.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import yaml

from src.core.logging import get_logger
from src.eval.tetramer_interface_analysis import (
    analyze_rescue_tetramer_impact,
    summarize_tetramer_analysis,
    CORE_DOMAIN_INTERFACE_RESIDUES,
)

logger = get_logger(__name__)


def run(args, configs: dict) -> int:
    """Run tetramer interface analysis on rescue candidates."""
    logger.info("Starting tetramer interface analysis...")

    # Paths
    paths_cfg = configs.get("paths", {})
    data_cfg = paths_cfg.get("data", {})
    processed_dir = Path(data_cfg.get("processed", "Data/processed"))
    raw_dir = Path(data_cfg.get("raw", "Data/raw"))
    reports_dir = Path(paths_cfg.get("reports_dir", "reports"))

    rescues_root = processed_dir / "rescues"
    tetramer_pdb = raw_dir / "experimental_pdbs" / "3KMD.pdb"

    if not tetramer_pdb.exists():
        logger.warning(f"Tetramer structure not found: {tetramer_pdb}")
        logger.info("Using pre-computed interface residues from 3KMD analysis")

    if not rescues_root.exists():
        logger.error(f"Rescues directory not found: {rescues_root}")
        return 1

    # Output directory
    tetramer_dir = reports_dir / "tetramer_interface"
    tetramer_dir.mkdir(parents=True, exist_ok=True)

    # Print interface residues
    print(f"\n{'='*80}")
    print(f"Tetramer Interface Residues (from 3KMD structure)")
    print(f"{'='*80}")
    print(f"Total interface positions: {len(CORE_DOMAIN_INTERFACE_RESIDUES)}")
    print(f"Positions: {CORE_DOMAIN_INTERFACE_RESIDUES}")
    print(f"\nKey interface regions:")
    print(f"  N-terminal (93-101): Dimer-dimer interface")
    print(f"  Zinc-binding (176-181, 242-244): Structural coordination")
    print(f"  H2 helix (198-202, 267): Dimer interface")
    print(f"  L3 loop (224-233): Conformational flexibility")

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

        logger.info(f"\nAnalyzing tetramer interface impact for {target}...")

        # Load rescue candidates
        pareto_df = pd.read_parquet(pareto_path)

        # Sort by Pareto status and quality
        pareto_df = pareto_df.sort_values(
            ["is_pareto", "n_rescue", "ddg_gain", "risk"],
            ascending=[False, True, True, True]
        )

        try:
            # Analyze tetramer interface impact
            tetramer_df = analyze_rescue_tetramer_impact(
                pareto_df.head(50),  # Analyze top 50
                tetramer_pdb=tetramer_pdb if tetramer_pdb.exists() else None,
                filter_high_risk=False,  # Don't filter yet, just annotate
            )

            # Save results
            output_csv = tetramer_dir / f"{target}_tetramer_interface.csv"
            tetramer_df.to_csv(output_csv, index=False)
            logger.info(f"Saved tetramer interface results to {output_csv}")

            # Generate summary
            summary = summarize_tetramer_analysis(
                tetramer_df,
                output_path=tetramer_dir / f"{target}_tetramer_summary.json"
            )

            # Print summary
            print(f"\n{'='*80}")
            print(f"Tetramer Interface Analysis: {target}")
            print(f"{'='*80}")
            print(f"  Rescues analyzed: {summary['n_rescues_analyzed']}")
            print(f"  Rescues at interface: {summary['n_at_interface']}")
            print(f"  Mean impact score: {summary['mean_impact_score']:.3f}")
            print(f"  Max impact score: {summary['max_impact_score']:.3f}")

            print(f"\n  Risk Distribution:")
            for risk_level, count in sorted(summary['risk_distribution'].items()):
                print(f"    {risk_level:12s}: {count}")

            # Identify dangerous rescues (high tetramer interface impact)
            dangerous = tetramer_df[
                tetramer_df["tetramer_risk"].isin(["critical", "high"])
            ]

            if not dangerous.empty:
                print(f"\n  ⚠️  WARNING: {len(dangerous)} rescues with high tetramer risk:")
                for idx, row in dangerous.head(5).iterrows():
                    print(f"    - {row['rescue_mutations']:20s} | "
                          f"Impact: {row['tetramer_impact_score']:.2f} | "
                          f"Risk: {row['tetramer_risk']}")
            else:
                print(f"\n  ✅ No high-risk tetramer interface conflicts detected")

            # Identify safe rescues
            safe = tetramer_df[
                tetramer_df["tetramer_risk"].isin(["safe", "low"])
            ]

            if not safe.empty:
                print(f"\n  ✅ Top 3 safe rescues (low tetramer impact):")
                for idx, row in safe.head(3).iterrows():
                    print(f"    - {row['rescue_mutations']:20s} | "
                          f"Impact: {row['tetramer_impact_score']:.2f} | "
                          f"Risk: {row['tetramer_risk']}")

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
        print("Consolidated Tetramer Interface Summary")
        print(f"{'='*80}")

        for result in all_results:
            target = result["target"]
            summary = result["summary"]
            print(f"\n{target}:")
            print(f"  Interface rescues: {summary['n_at_interface']}")
            print(f"  High-risk rescues: {result['n_dangerous']}")
            print(f"  Safe rescues: {result['n_safe']}")
            print(f"  Mean impact score: {summary['mean_impact_score']:.3f}")

        # Save consolidated CSV
        all_tetramer = []
        for target_dir in sorted(rescues_root.iterdir()):
            if not target_dir.is_dir():
                continue
            target = target_dir.name
            csv_path = tetramer_dir / f"{target}_tetramer_interface.csv"
            if csv_path.exists():
                df = pd.read_csv(csv_path)
                df["target"] = target
                all_tetramer.append(df)

        if all_tetramer:
            consolidated_df = pd.concat(all_tetramer, ignore_index=True)
            consolidated_csv = tetramer_dir / "all_targets_tetramer_interface.csv"
            consolidated_df.to_csv(consolidated_csv, index=False)
            logger.info(f"Saved consolidated results to {consolidated_csv}")

    print(f"\n{'='*80}")
    print(f"Tetramer interface analysis complete!")
    print(f"Results saved to: {tetramer_dir}")
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
