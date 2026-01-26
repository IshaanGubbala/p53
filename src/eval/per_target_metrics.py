"""Per-target (per-position) validation metrics for hotspot separation analysis."""
from __future__ import annotations

import logging
from typing import Optional

import numpy as np
import pandas as pd

from .variant_separation import compute_auc, bootstrap_ci


logger = logging.getLogger(__name__)


def compute_per_target_auc(
    df: pd.DataFrame,
    min_variants: int = 5,
) -> pd.DataFrame:
    """Compute AUC per position (target hotspot).

    Args:
        df: DataFrame with columns: pos, ddg, label
        min_variants: Minimum variants per position to compute AUC

    Returns:
        DataFrame with columns: pos, n_variants, n_benign, n_pathogenic, auc
    """
    if not {"pos", "ddg", "label"}.issubset(df.columns):
        raise ValueError("DataFrame missing required columns: pos, ddg, label")

    results: list[dict] = []

    for pos, group in df.groupby("pos"):
        n_variants = len(group)
        if n_variants < min_variants:
            continue

        n_benign = (group["label"] == 0).sum()
        n_pathogenic = (group["label"] == 1).sum()

        # Need both classes to compute AUC
        if n_benign == 0 or n_pathogenic == 0:
            continue

        try:
            auc = compute_auc(group["label"], group["ddg"])
            results.append({
                "pos": int(pos),
                "n_variants": int(n_variants),
                "n_benign": int(n_benign),
                "n_pathogenic": int(n_pathogenic),
                "auc": float(auc),
            })
        except Exception as exc:
            logger.warning(f"Failed to compute AUC for position {pos}: {exc}")
            continue

    if not results:
        logger.warning("No positions with sufficient variants for per-target AUC")
        return pd.DataFrame(columns=["pos", "n_variants", "n_benign", "n_pathogenic", "auc"])

    return pd.DataFrame(results).sort_values("pos").reset_index(drop=True)


def per_target_bootstrap_ci(
    df: pd.DataFrame,
    min_variants: int = 5,
    n_bootstrap: int = 2000,
    confidence_level: float = 0.95,
    random_seed: int = 42,
) -> pd.DataFrame:
    """Compute per-position AUC with bootstrap confidence intervals.

    Args:
        df: DataFrame with columns: pos, ddg, label
        min_variants: Minimum variants per position
        n_bootstrap: Number of bootstrap iterations
        confidence_level: Confidence level (e.g., 0.95 for 95% CI)
        random_seed: Random seed for reproducibility

    Returns:
        DataFrame with columns: pos, n_variants, auc, ci_lower, ci_upper
    """
    if not {"pos", "ddg", "label"}.issubset(df.columns):
        raise ValueError("DataFrame missing required columns: pos, ddg, label")

    # Compute quantiles from confidence level
    alpha = 1.0 - confidence_level
    lower_quantile = alpha / 2
    upper_quantile = 1.0 - alpha / 2

    results: list[dict] = []

    for pos, group in df.groupby("pos"):
        n_variants = len(group)
        if n_variants < min_variants:
            continue

        n_benign = (group["label"] == 0).sum()
        n_pathogenic = (group["label"] == 1).sum()

        if n_benign == 0 or n_pathogenic == 0:
            continue

        try:
            # Compute point estimate
            auc = compute_auc(group["label"], group["ddg"])

            # Bootstrap CI
            ci_lower, ci_upper = bootstrap_ci(
                group["label"],
                group["ddg"],
                n=n_bootstrap,
                seed=random_seed,
            )

            results.append({
                "pos": int(pos),
                "n_variants": int(n_variants),
                "n_benign": int(n_benign),
                "n_pathogenic": int(n_pathogenic),
                "auc": float(auc),
                "ci_lower": float(ci_lower),
                "ci_upper": float(ci_upper),
                "ci_width": float(ci_upper - ci_lower),
            })
        except Exception as exc:
            logger.warning(f"Failed to compute bootstrap CI for position {pos}: {exc}")
            continue

    if not results:
        logger.warning("No positions with sufficient variants for bootstrap CI")
        return pd.DataFrame(columns=[
            "pos", "n_variants", "n_benign", "n_pathogenic",
            "auc", "ci_lower", "ci_upper", "ci_width"
        ])

    return pd.DataFrame(results).sort_values("pos").reset_index(drop=True)


def filter_hotspot_positions(
    df: pd.DataFrame,
    hotspot_positions: list[int],
) -> pd.DataFrame:
    """Filter DataFrame to only include hotspot positions.

    Args:
        df: DataFrame with 'pos' column
        hotspot_positions: List of positions to include

    Returns:
        Filtered DataFrame
    """
    if "pos" not in df.columns:
        raise ValueError("DataFrame missing 'pos' column")

    return df[df["pos"].isin(hotspot_positions)].copy()


def compare_train_test_per_target(
    train_df: pd.DataFrame,
    test_df: pd.DataFrame,
    min_variants: int = 5,
    n_bootstrap: int = 2000,
    random_seed: int = 42,
) -> pd.DataFrame:
    """Compare per-target AUC between train and test sets.

    Args:
        train_df: Training set with pos, ddg, label
        test_df: Test set with pos, ddg, label
        min_variants: Minimum variants per position
        n_bootstrap: Bootstrap iterations
        random_seed: Random seed

    Returns:
        DataFrame with columns: pos, train_auc, train_ci_lower, train_ci_upper,
                                 test_auc, test_ci_lower, test_ci_upper
    """
    train_results = per_target_bootstrap_ci(
        train_df,
        min_variants=min_variants,
        n_bootstrap=n_bootstrap,
        random_seed=random_seed,
    )

    test_results = per_target_bootstrap_ci(
        test_df,
        min_variants=min_variants,
        n_bootstrap=n_bootstrap,
        random_seed=random_seed + 1,  # Different seed for test
    )

    # Merge on position
    if train_results.empty or test_results.empty:
        return pd.DataFrame()

    merged = train_results.merge(
        test_results[["pos", "n_variants", "auc", "ci_lower", "ci_upper"]],
        on="pos",
        how="outer",
        suffixes=("_train", "_test"),
    )

    # Reorder columns
    column_order = [
        "pos",
        "n_variants_train", "auc_train", "ci_lower_train", "ci_upper_train",
        "n_variants_test", "auc_test", "ci_lower_test", "ci_upper_test",
    ]
    available_columns = [col for col in column_order if col in merged.columns]

    return merged[available_columns].sort_values("pos").reset_index(drop=True)


def identify_challenging_positions(
    per_target_df: pd.DataFrame,
    auc_threshold: float = 0.7,
    ci_width_threshold: float = 0.15,
) -> pd.DataFrame:
    """Identify positions with poor separation or high uncertainty.

    Args:
        per_target_df: DataFrame from per_target_bootstrap_ci()
        auc_threshold: AUC below this is considered challenging
        ci_width_threshold: CI width above this is considered uncertain

    Returns:
        DataFrame of challenging positions with reasons
    """
    if per_target_df.empty:
        return pd.DataFrame(columns=["pos", "auc", "ci_width", "reason"])

    challenging: list[dict] = []

    for _, row in per_target_df.iterrows():
        reasons = []
        if row["auc"] < auc_threshold:
            reasons.append(f"low_auc_{row['auc']:.3f}")
        if row["ci_width"] > ci_width_threshold:
            reasons.append(f"wide_ci_{row['ci_width']:.3f}")

        if reasons:
            challenging.append({
                "pos": int(row["pos"]),
                "auc": float(row["auc"]),
                "ci_width": float(row["ci_width"]),
                "n_variants": int(row["n_variants"]),
                "reason": ",".join(reasons),
            })

    if not challenging:
        return pd.DataFrame(columns=["pos", "auc", "ci_width", "n_variants", "reason"])

    return pd.DataFrame(challenging).sort_values("auc").reset_index(drop=True)
