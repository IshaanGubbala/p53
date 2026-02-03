"""
Run druggability analysis on top rescue candidates.

Alternative to fpocket that works with existing SASA and structural data.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

from src.core.logging import get_logger
from src.eval.druggability_analysis import (
    analyze_rescue_druggability,
    summarize_druggability,
)
from src.viz.plots_druggability import (
    plot_druggability_distribution,
    plot_druggability_comparison,
    plot_top_druggable_rescues,
)

logger = get_logger(__name__)


def run(args, configs: dict) -> int:
    """Run druggability analysis on rescue candidates."""
    logger.info("Starting druggability analysis...")

    # Paths
    paths_cfg = configs.get("paths", {})
    data_cfg = paths_cfg.get("data", {})
    processed_dir = Path(data_cfg.get("processed", "Data/processed"))
    interim_dir = Path(data_cfg.get("interim", "Data/interim"))
    reports_dir = Path(paths_cfg.get("reports_dir", "reports"))

    rescues_root = processed_dir / "rescues"
    base_pdb = interim_dir / "structure" / "P04637_core_94_312.pdb"
    sasa_path = interim_dir / "structure" / "P04637_core_94_312_sasa.json"

    # EvoEF2 config for building rescued structures
    # Note: scoring.yaml has top-level keys, not nested under "scoring"
    evoef2_cfg = configs.get("evoef2", {})
    work_dir = processed_dir / "cache" / "druggability_structures"

    # Load functional sites for filtering
    functional_sites = configs.get("functional_sites", {})
    if functional_sites:
        total_functional = sum(len(sites) for sites in functional_sites.values())
        logger.info(f"Loaded {total_functional} functional sites across {len(functional_sites)} categories")
        logger.info(f"Will filter rescues that mutate: {', '.join(functional_sites.keys())}")
    else:
        logger.warning("No functional sites loaded - proceeding without filtering")

    if not base_pdb.exists():
        logger.error(f"Base PDB not found: {base_pdb}")
        return 1

    if not rescues_root.exists():
        logger.error(f"Rescues directory not found: {rescues_root}")
        return 1

    # Output directory
    druggability_dir = reports_dir / "druggability"
    druggability_dir.mkdir(parents=True, exist_ok=True)

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

        logger.info(f"Analyzing druggability for {target}...")

        # Load rescue candidates
        pareto_df = pd.read_parquet(pareto_path)

        # Sort by Pareto status and quality
        pareto_df = pareto_df.sort_values(
            ["is_pareto", "n_rescue", "ddg_gain", "risk"],
            ascending=[False, True, True, True]
        )

        try:
            # Analyze druggability by BUILDING rescued structures
            # This is computationally expensive but necessary for accurate analysis
            logger.info(f"Building rescued structures for {target} (this may take a few minutes)...")

            druggability_df = analyze_rescue_druggability(
                pareto_df,
                target=target,
                base_pdb=base_pdb,  # Wild-type structure (template for building)
                sasa_path=sasa_path,
                top_n=50,  # Analyze more candidates since we're filtering
                evoef2_cfg=evoef2_cfg,  # For BuildMutant
                work_dir=work_dir,  # Cache directory
                functional_sites=functional_sites,  # Filter critical sites
                filter_functional=True,  # Enable filtering
            )

            # Save results
            output_csv = druggability_dir / f"{target}_druggability.csv"
            druggability_df.to_csv(output_csv, index=False)
            logger.info(f"Saved druggability results to {output_csv}")

            # Generate visualizations
            try:
                plot_druggability_distribution(
                    druggability_df,
                    druggability_dir / f"{target}_druggability_dist.png",
                    target=target,
                )
                plot_top_druggable_rescues(
                    druggability_df,
                    druggability_dir / f"{target}_top_druggable.png",
                    target=target,
                    top_n=10,
                )
            except Exception as e:
                logger.warning(f"Failed to generate plots for {target}: {e}")

            # Generate summary
            summary = summarize_druggability(
                druggability_df,
                output_path=druggability_dir / f"{target}_druggability_summary.json"
            )

            # Count functional site conflicts (if any remain)
            n_conflicts = druggability_df.get("functional_site_conflict", pd.Series([False])).sum()

            # Print summary
            print(f"\n{'='*80}")
            print(f"Druggability Analysis: {target}")
            print(f"{'='*80}")
            print(f"  Rescues analyzed: {summary['n_rescues_analyzed']}")
            if functional_sites:
                print(f"  ⚠️  Filtered rescues (functional sites): many (see logs)")
                print(f"  ℹ️  Remaining results are 'safe' - no critical site mutations")
            print(f"  High druggability (≥0.6): {summary['high_druggability']}")
            print(f"  Medium druggability (0.4-0.6): {summary['medium_druggability']}")
            print(f"  Low druggability (<0.4): {summary['low_druggability']}")
            print(f"  Mean score: {summary['mean_score']}")
            print(f"  Max score: {summary['max_score']}")

            print(f"\n  Top 3 druggable rescues:")
            for i, rescue in enumerate(summary['top_druggable_rescues'], 1):
                print(f"    {i}. {rescue['rescue_mutations']:20s} | "
                      f"Score: {rescue['druggability_score']:.3f} | "
                      f"{rescue.get('pocket_flag', 'N/A')}")

            all_results.append({
                "target": target,
                "summary": summary,
                "top_rescues": druggability_df.head(5),
            })

        except Exception as e:
            logger.error(f"Failed to analyze {target}: {e}", exc_info=True)
            continue

    # Generate consolidated report
    if all_results:
        print(f"\n{'='*80}")
        print("Consolidated Druggability Summary")
        print(f"{'='*80}")

        for result in all_results:
            target = result["target"]
            summary = result["summary"]
            print(f"\n{target}:")
            print(f"  High/Medium/Low: {summary['high_druggability']}/"
                  f"{summary['medium_druggability']}/{summary['low_druggability']}")
            print(f"  Best rescue: {summary['top_druggable_rescues'][0]['rescue_mutations']} "
                  f"(score: {summary['top_druggable_rescues'][0]['druggability_score']})")

        # Save consolidated CSV
        all_druggable = []
        for result in all_results:
            top_df = result["top_rescues"].copy()
            top_df["target"] = result["target"]
            all_druggable.append(top_df)

        if all_druggable:
            consolidated_df = pd.concat(all_druggable, ignore_index=True)
            consolidated_csv = druggability_dir / "all_targets_druggability.csv"
            consolidated_df.to_csv(consolidated_csv, index=False)
            logger.info(f"Saved consolidated results to {consolidated_csv}")

            # Generate comparison plot
            try:
                plot_druggability_comparison(
                    consolidated_df,
                    druggability_dir / "druggability_comparison.png",
                )
            except Exception as e:
                logger.warning(f"Failed to generate comparison plot: {e}")

    print(f"\n{'='*80}")
    print(f"Druggability analysis complete!")
    print(f"Results saved to: {druggability_dir}")
    print(f"{'='*80}\n")

    return 0


if __name__ == "__main__":
    import yaml
    from pathlib import Path

    # Load configs
    config_files = [
        "configs/p53.yaml",
        "configs/paths.yaml",
        "configs/scoring.yaml",  # Need EvoEF2 config
    ]

    configs = {}
    for cfg_file in config_files:
        cfg_path = Path(cfg_file)
        if cfg_path.exists():
            with open(cfg_path) as f:
                cfg_data = yaml.safe_load(f) or {}
                configs.update(cfg_data)

    sys.exit(run(None, configs))
