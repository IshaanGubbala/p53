"""
Risk weight sensitivity analysis for rescue mutations.

This module tests the robustness of rescue rankings under different
risk weighting schemes to ensure results aren't artifacts of parameter tuning.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from src.core.logging import get_logger

logger = get_logger(__name__)


# Define alternative weighting profiles
WEIGHT_PROFILES = {
    "balanced": {
        "name": "Balanced (Current)",
        "description": "Balance all factors with slight bias toward functional sites",
        "weights": {
            "functional": 0.40,
            "conservation": 0.20,
            "burial": 0.15,
            "msa_conservation": 0.25,
        },
    },
    "safety_first": {
        "name": "Safety-First",
        "description": "Prioritize not breaking existing function above all else",
        "weights": {
            "functional": 0.70,
            "conservation": 0.20,
            "burial": 0.05,
            "msa_conservation": 0.05,
        },
    },
    "evolution_pure": {
        "name": "Evolution-Pure",
        "description": "Trust evolution - highly conserved positions are conserved for a reason",
        "weights": {
            "functional": 0.10,
            "conservation": 0.30,
            "burial": 0.10,
            "msa_conservation": 0.50,
        },
    },
    "structure_first": {
        "name": "Structure-First",
        "description": "Prioritize restoring protein stability/folding",
        "weights": {
            "functional": 0.20,
            "conservation": 0.10,
            "burial": 0.50,
            "msa_conservation": 0.20,
        },
    },
}


def calculate_risk_with_weights(
    df: pd.DataFrame,
    weights: dict[str, float],
    risk_columns: dict[str, str] | None = None,
) -> pd.Series:
    """
    Calculate risk score using specified weights.

    Args:
        df: DataFrame with risk component columns
        weights: Dict of component -> weight
        risk_columns: Optional dict mapping weight key to actual column name

    Returns:
        Series with calculated risk scores
    """
    if risk_columns is None:
        risk_columns = {
            "functional": "functional_risk",
            "conservation": "conservation_risk",
            "burial": "burial_risk",
            "msa_conservation": "msa_conservation_risk",
        }

    risk_score = pd.Series(0.0, index=df.index)

    for component, weight in weights.items():
        col_name = risk_columns.get(component, f"{component}_risk")
        if col_name in df.columns:
            risk_score += weight * df[col_name]
        else:
            logger.warning(f"Column {col_name} not found, skipping component {component}")

    return risk_score


def sensitivity_analysis(
    variants_df: pd.DataFrame,
    weight_profiles: dict[str, dict] | None = None,
    top_n: int = 10,
    risk_columns: dict[str, str] | None = None,
) -> dict[str, Any]:
    """
    Perform risk weight sensitivity analysis.

    Args:
        variants_df: DataFrame with all risk components
        weight_profiles: Dict of profile_name -> weight_dict (uses WEIGHT_PROFILES if None)
        top_n: Number of top candidates to track
        risk_columns: Optional dict mapping weight key to actual column name

    Returns:
        Dict with analysis results
    """
    if weight_profiles is None:
        weight_profiles = WEIGHT_PROFILES

    logger.info(f"Running sensitivity analysis with {len(weight_profiles)} weighting profiles")

    # Store results for each profile
    profile_results = {}
    all_rankings = {}

    for profile_id, profile_config in weight_profiles.items():
        logger.info(f"Analyzing profile: {profile_config['name']}")

        weights = profile_config["weights"]

        # Calculate risk scores with this profile
        risk_scores = calculate_risk_with_weights(variants_df, weights, risk_columns)

        # Rank candidates (lower risk = better rank)
        variant_ids = variants_df.get("rescue_mutations", variants_df.index)
        rankings = risk_scores.rank(method="min").astype(int)

        # Get top N
        top_indices = risk_scores.nsmallest(top_n).index
        top_variants = variant_ids.iloc[top_indices].tolist()

        profile_results[profile_id] = {
            "name": profile_config["name"],
            "description": profile_config["description"],
            "weights": weights,
            "top_candidates": top_variants,
            "top_risk_scores": risk_scores.iloc[top_indices].tolist(),
        }

        # Store rankings for all variants
        for idx, variant_id in zip(variant_ids.index, variant_ids):
            if variant_id not in all_rankings:
                all_rankings[variant_id] = {}
            all_rankings[variant_id][profile_id] = int(rankings.loc[idx])

    # Compute robustness metrics
    logger.info("Computing robustness metrics...")

    robustness_scores = {}
    for variant_id, rankings_dict in all_rankings.items():
        # Count how many profiles put this variant in top N
        in_top_n_count = sum(1 for rank in rankings_dict.values() if rank <= top_n)

        # Calculate rank variance (lower = more consistent)
        rank_variance = float(np.var(list(rankings_dict.values())))

        # Calculate median rank
        median_rank = float(np.median(list(rankings_dict.values())))

        # Robustness score: fraction of profiles where variant is in top N
        robustness_score = in_top_n_count / len(weight_profiles)

        robustness_scores[variant_id] = {
            "in_top_n_count": in_top_n_count,
            "rank_variance": rank_variance,
            "median_rank": median_rank,
            "robustness_score": robustness_score,
            "rankings": rankings_dict,
        }

    # Identify highly robust candidates (in top N for all profiles)
    highly_robust = [
        variant_id
        for variant_id, metrics in robustness_scores.items()
        if metrics["robustness_score"] == 1.0
    ]

    # Identify robust candidates (in top N for ≥3/4 profiles)
    robust = [
        variant_id
        for variant_id, metrics in robustness_scores.items()
        if metrics["robustness_score"] >= 0.75 and metrics["robustness_score"] < 1.0
    ]

    # Identify moderately robust (in top N for 2/4 profiles)
    moderately_robust = [
        variant_id
        for variant_id, metrics in robustness_scores.items()
        if metrics["robustness_score"] >= 0.50 and metrics["robustness_score"] < 0.75
    ]

    return {
        "profile_results": profile_results,
        "robustness_scores": robustness_scores,
        "highly_robust": highly_robust,
        "robust": robust,
        "moderately_robust": moderately_robust,
        "n_profiles": len(weight_profiles),
        "top_n": top_n,
    }


def generate_sensitivity_report(
    sensitivity_results: dict[str, Any],
    output_dir: Path,
    target: str,
) -> Path:
    """
    Generate sensitivity analysis report with rankings table.

    Args:
        sensitivity_results: Results from sensitivity_analysis()
        output_dir: Directory to save report
        target: Target mutation name (e.g., "R175H")

    Returns:
        Path to generated report
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Create rankings DataFrame
    robustness_scores = sensitivity_results["robustness_scores"]
    profile_results = sensitivity_results["profile_results"]

    rows = []
    for variant_id, metrics in robustness_scores.items():
        row = {
            "rescue_mutations": variant_id,
            "robustness_score": metrics["robustness_score"],
            "median_rank": metrics["median_rank"],
            "rank_variance": metrics["rank_variance"],
        }

        # Add rankings for each profile
        for profile_id in profile_results.keys():
            row[f"rank_{profile_id}"] = metrics["rankings"].get(profile_id, None)

        rows.append(row)

    rankings_df = pd.DataFrame(rows)

    # Sort by robustness score (descending) then median rank (ascending)
    rankings_df = rankings_df.sort_values(
        ["robustness_score", "median_rank"],
        ascending=[False, True]
    )

    # Save CSV
    csv_path = output_dir / f"{target}_sensitivity_rankings.csv"
    rankings_df.to_csv(csv_path, index=False)
    logger.info(f"Saved sensitivity rankings to {csv_path}")

    # Save JSON summary
    json_path = output_dir / f"{target}_sensitivity_summary.json"

    summary = {
        "target": target,
        "n_profiles": sensitivity_results["n_profiles"],
        "top_n": sensitivity_results["top_n"],
        "highly_robust": sensitivity_results["highly_robust"],
        "robust": sensitivity_results["robust"],
        "moderately_robust": sensitivity_results["moderately_robust"],
        "profile_descriptions": {
            pid: {
                "name": pdata["name"],
                "description": pdata["description"],
                "weights": pdata["weights"],
                "top_3": pdata["top_candidates"][:3],
            }
            for pid, pdata in profile_results.items()
        },
    }

    with open(json_path, "w") as f:
        json.dump(summary, f, indent=2)

    logger.info(f"Saved sensitivity summary to {json_path}")

    return csv_path


def classify_robustness(robustness_score: float) -> str:
    """
    Classify robustness score into categories.

    Args:
        robustness_score: Fraction of profiles where variant is in top N (0-1)

    Returns:
        Classification string
    """
    if robustness_score == 1.0:
        return "HIGHLY ROBUST"
    elif robustness_score >= 0.75:
        return "ROBUST"
    elif robustness_score >= 0.50:
        return "MODERATELY ROBUST"
    else:
        return "SENSITIVE"
