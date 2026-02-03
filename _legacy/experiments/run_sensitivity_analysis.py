"""
Run risk weight sensitivity analysis on rescue candidates.

This analysis tests robustness of rescue rankings under different
risk weighting schemes.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import yaml

from src.core.logging import get_logger
from src.eval.sensitivity_analysis import (
    sensitivity_analysis,
    generate_sensitivity_report,
    classify_robustness,
    WEIGHT_PROFILES,
)

logger = get_logger(__name__)


def run(args, configs: dict) -> int:
    """Run sensitivity analysis on rescue candidates."""
    logger.info("Starting risk weight sensitivity analysis...")

    # Paths
    paths_cfg = configs.get("paths", {})
    data_cfg = paths_cfg.get("data", {})
    processed_dir = Path(data_cfg.get("processed", "Data/processed"))
    reports_dir = Path(paths_cfg.get("reports_dir", "reports"))

    rescues_root = processed_dir / "rescues"

    if not rescues_root.exists():
        logger.error(f"Rescues directory not found: {rescues_root}")
        return 1

    # Output directory
    sensitivity_dir = reports_dir / "sensitivity"
    sensitivity_dir.mkdir(parents=True, exist_ok=True)

    # Print weighting profiles
    print(f"\n{'='*80}")
    print(f"Risk Weighting Profiles")
    print(f"{'='*80}")
    for profile_id, profile_config in WEIGHT_PROFILES.items():
        print(f"\n{profile_config['name']}:")
        print(f"  {profile_config['description']}")
        print(f"  Weights:")
        for component, weight in profile_config['weights'].items():
            print(f"    {component:20s}: {weight:.2f}")

    # Process each target
    all_results = []

    for target_dir in sorted(rescues_root.iterdir()):
        if not target_dir.is_dir():
            continue

        target = target_dir.name
        candidates_path = target_dir / "candidates.parquet"

        if not candidates_path.exists():
            logger.warning(f"No candidates file for {target}, skipping")
            continue

        logger.info(f"\nRunning sensitivity analysis for {target}...")

        # Load rescue candidates
        candidates_df = pd.read_parquet(candidates_path)

        # Expand risk_components JSON column into separate columns
        if "risk_components" in candidates_df.columns:
            import json
            risk_components = candidates_df["risk_components"].apply(
                lambda x: json.loads(x) if isinstance(x, str) else x
            )
            risk_df = pd.DataFrame(risk_components.tolist())

            # Rename to match expected column names
            rename_map = {
                "functional": "functional_risk",
                "conservation": "conservation_risk",
                "burial": "burial_risk",
                "msa_conservation": "msa_conservation_risk",
            }
            risk_df = risk_df.rename(columns=rename_map)

            # Add to candidates_df
            for col in risk_df.columns:
                candidates_df[col] = risk_df[col]

        # Check required columns
        required_cols = ["functional_risk", "conservation_risk", "burial_risk"]
        missing_cols = [col for col in required_cols if col not in candidates_df.columns]

        if missing_cols:
            logger.warning(f"Missing risk columns for {target}: {missing_cols}, skipping")
            continue

        # Add MSA conservation risk if available (otherwise use 0)
        if "msa_conservation_risk" not in candidates_df.columns:
            logger.warning(f"No MSA conservation risk for {target}, using 0")
            candidates_df["msa_conservation_risk"] = 0.0

        try:
            # Run sensitivity analysis
            sensitivity_results = sensitivity_analysis(
                candidates_df,
                weight_profiles=WEIGHT_PROFILES,
                top_n=10,
            )

            # Generate report
            generate_sensitivity_report(
                sensitivity_results,
                sensitivity_dir,
                target,
            )

            # Print summary
            print(f"\n{'='*80}")
            print(f"Sensitivity Analysis: {target}")
            print(f"{'='*80}")
            print(f"  Candidates analyzed: {len(candidates_df)}")
            print(f"  Weighting profiles: {sensitivity_results['n_profiles']}")
            print(f"  Top N tracked: {sensitivity_results['top_n']}")

            print(f"\n  Robustness Summary:")
            print(f"    Highly Robust (top 10 in ALL profiles): {len(sensitivity_results['highly_robust'])}")
            print(f"    Robust (top 10 in ≥3/4 profiles): {len(sensitivity_results['robust'])}")
            print(f"    Moderately Robust (top 10 in 2/4): {len(sensitivity_results['moderately_robust'])}")

            # Print highly robust candidates
            if sensitivity_results['highly_robust']:
                print(f"\n  🌟 HIGHLY ROBUST Rescues (top 10 in ALL profiles):")
                for variant_id in sensitivity_results['highly_robust'][:5]:
                    metrics = sensitivity_results['robustness_scores'][variant_id]
                    print(f"    - {variant_id:25s} | "
                          f"Median rank: {metrics['median_rank']:.0f}")
            else:
                print(f"\n  ⚠️  No highly robust rescues found (none in top 10 for all profiles)")

            # Print top candidates for each profile
            print(f"\n  Top 3 Candidates by Profile:")
            for profile_id, profile_data in sensitivity_results['profile_results'].items():
                print(f"\n    {profile_data['name']}:")
                for i, variant_id in enumerate(profile_data['top_candidates'][:3], 1):
                    risk_score = profile_data['top_risk_scores'][i-1]
                    print(f"      {i}. {variant_id:25s} | Risk: {risk_score:.3f}")

            all_results.append({
                "target": target,
                "highly_robust": sensitivity_results['highly_robust'],
                "robust": sensitivity_results['robust'],
                "n_candidates": len(candidates_df),
            })

        except Exception as e:
            logger.error(f"Failed to analyze {target}: {e}", exc_info=True)
            continue

    # Generate consolidated report
    if all_results:
        print(f"\n{'='*80}")
        print("Consolidated Sensitivity Summary")
        print(f"{'='*80}")

        for result in all_results:
            target = result["target"]
            print(f"\n{target}:")
            print(f"  Total candidates: {result['n_candidates']}")
            print(f"  Highly robust: {len(result['highly_robust'])}")
            print(f"  Robust: {len(result['robust'])}")

            if result['highly_robust']:
                print(f"  Top highly robust: {', '.join(result['highly_robust'][:3])}")

        # Identify universal winners (highly robust across ALL targets)
        all_highly_robust = {}
        for result in all_results:
            for variant_id in result['highly_robust']:
                if variant_id not in all_highly_robust:
                    all_highly_robust[variant_id] = 0
                all_highly_robust[variant_id] += 1

        universal_winners = [
            variant_id
            for variant_id, count in all_highly_robust.items()
            if count >= len(all_results)  # Appears in all targets
        ]

        if universal_winners:
            print(f"\n{'='*80}")
            print(f"🏆 UNIVERSAL WINNERS (highly robust across ALL targets):")
            print(f"{'='*80}")
            for variant_id in universal_winners:
                print(f"  - {variant_id}")
        else:
            print(f"\n(No universal winners - rescues are target-specific)")

    print(f"\n{'='*80}")
    print(f"Sensitivity analysis complete!")
    print(f"Results saved to: {sensitivity_dir}")
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
