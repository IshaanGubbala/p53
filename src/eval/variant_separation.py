from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


def _make_mutation_key(df: pd.DataFrame) -> pd.Series:
    return df["ref"].astype(str).str.upper() + df["pos"].astype(int).astype(str) + df["alt"].astype(str).str.upper()


def load_labelled_scores(scores_path: Path, labels_dir: Path) -> pd.DataFrame:
    scores = pd.read_parquet(scores_path)
    if not {"pos", "ref", "alt", "ddg"}.issubset(scores.columns):
        raise RuntimeError("Variant scores missing required columns: pos/ref/alt/ddg")

    scores = scores.copy()
    scores["pos"] = scores["pos"].astype(int)
    scores["ref"] = scores["ref"].astype(str).str.upper()
    scores["alt"] = scores["alt"].astype(str).str.upper()
    scores["mutation"] = _make_mutation_key(scores)
    scores = scores.groupby(["pos", "ref", "alt"], as_index=False)["ddg"].mean()

    benign = pd.read_parquet(labels_dir / "benign.parquet")
    pathogenic = pd.read_parquet(labels_dir / "pathogenic.parquet")
    if benign.empty or pathogenic.empty:
        raise RuntimeError("Label sets are empty; rerun build-dataset with ClinVar enabled")

    benign = benign.copy()
    pathogenic = pathogenic.copy()
    benign["label"] = 0
    pathogenic["label"] = 1
    labels = pd.concat([benign, pathogenic], ignore_index=True)
    labels["pos"] = labels["pos"].astype(int)
    labels["ref"] = labels["ref"].astype(str).str.upper()
    labels["alt"] = labels["alt"].astype(str).str.upper()

    grouped = labels.groupby(["pos", "ref", "alt"])["label"].nunique()
    conflict_keys = grouped[grouped > 1].index
    if len(conflict_keys) > 0:
        labels = labels.set_index(["pos", "ref", "alt"])
        labels = labels.loc[~labels.index.isin(conflict_keys)].reset_index()

    labels = labels.drop_duplicates(subset=["pos", "ref", "alt"], keep="first")

    merged = scores.merge(labels[["pos", "ref", "alt", "label"]], on=["pos", "ref", "alt"], how="inner")
    merged["mutation"] = _make_mutation_key(merged)
    return merged


def compute_auc(labels: Iterable[int], scores: Iterable[float]) -> float:
    labels_arr = np.asarray(list(labels))
    scores_arr = np.asarray(list(scores))
    if labels_arr.ndim != 1 or scores_arr.ndim != 1:
        raise ValueError("labels and scores must be 1D")
    if len(labels_arr) != len(scores_arr):
        raise ValueError("labels and scores must have same length")

    pos_mask = labels_arr == 1
    neg_mask = labels_arr == 0
    n_pos = int(pos_mask.sum())
    n_neg = int(neg_mask.sum())
    if n_pos == 0 or n_neg == 0:
        raise ValueError("Both positive and negative labels are required for AUC")

    ranks = pd.Series(scores_arr).rank(method="average").to_numpy()
    sum_ranks_pos = ranks[pos_mask].sum()
    auc = (sum_ranks_pos - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)
    return float(auc)


def bootstrap_ci(
    labels: Iterable[int],
    scores: Iterable[float],
    n: int = 2000,
    seed: int = 1337,
) -> tuple[float, float]:
    """
    Compute bootstrap confidence interval for AUC.

    NOTE: This uses simple bootstrap without stratification.
    For imbalanced datasets, use bootstrap_ci_stratified() instead.
    """
    labels_arr = np.asarray(list(labels))
    scores_arr = np.asarray(list(scores))
    rng = np.random.default_rng(seed)
    n_samples = len(labels_arr)
    aucs: list[float] = []
    for _ in range(n):
        idx = rng.integers(0, n_samples, size=n_samples)
        aucs.append(compute_auc(labels_arr[idx], scores_arr[idx]))
    lo, hi = np.quantile(aucs, [0.025, 0.975])
    return float(lo), float(hi)


def bootstrap_ci_stratified(
    labels: Iterable[int],
    scores: Iterable[float],
    n: int = 10000,
    seed: int = 1337,
    confidence: float = 0.95,
) -> tuple[float, float, float]:
    """
    Compute STRATIFIED bootstrap confidence interval for AUC.

    This is the CORRECT method for imbalanced datasets because:
    - Simple bootstrap may create samples with no positives or negatives
    - Stratified bootstrap maintains class ratios in each sample
    - More accurate confidence intervals, especially for small datasets

    Args:
        labels: Binary labels (0/1)
        scores: Prediction scores (higher = more likely positive)
        n: Number of bootstrap iterations (10000 recommended)
        seed: Random seed
        confidence: Confidence level (default 0.95 for 95% CI)

    Returns:
        Tuple of (lower_bound, upper_bound, point_estimate)
    """
    labels_arr = np.asarray(list(labels))
    scores_arr = np.asarray(list(scores))
    rng = np.random.default_rng(seed)

    # Get indices for each class
    pos_idx = np.where(labels_arr == 1)[0]
    neg_idx = np.where(labels_arr == 0)[0]

    n_pos = len(pos_idx)
    n_neg = len(neg_idx)

    if n_pos == 0 or n_neg == 0:
        raise ValueError("Both positive and negative labels are required")

    aucs: list[float] = []

    for _ in range(n):
        # Stratified sampling: sample with replacement within each class
        boot_pos = rng.choice(pos_idx, size=n_pos, replace=True)
        boot_neg = rng.choice(neg_idx, size=n_neg, replace=True)

        # Combine
        boot_idx = np.concatenate([boot_pos, boot_neg])

        try:
            auc = compute_auc(labels_arr[boot_idx], scores_arr[boot_idx])
            aucs.append(auc)
        except ValueError:
            # Skip invalid samples (shouldn't happen with stratified)
            continue

    if len(aucs) < n * 0.9:
        raise RuntimeError(f"Too many invalid bootstrap samples: {n - len(aucs)}/{n}")

    # Compute percentiles
    alpha = (1 - confidence) / 2
    lo, hi = np.quantile(aucs, [alpha, 1 - alpha])
    point_estimate = compute_auc(labels_arr, scores_arr)

    return float(lo), float(hi), float(point_estimate)


def compute_auprc(labels: Iterable[int], scores: Iterable[float]) -> float:
    """
    Compute Area Under Precision-Recall Curve.

    AUPRC is more informative than AUC for imbalanced datasets because:
    - AUC can be high even when many positives are missed
    - AUPRC directly measures precision at each recall level
    - More sensitive to performance on the minority class

    Args:
        labels: Binary labels (0/1)
        scores: Prediction scores (higher = more likely positive)

    Returns:
        AUPRC score (0-1, higher is better)
    """
    labels_arr = np.asarray(list(labels))
    scores_arr = np.asarray(list(scores))

    # Sort by descending scores
    desc_order = np.argsort(-scores_arr)
    labels_sorted = labels_arr[desc_order]

    # Compute precision at each threshold
    n_pos = labels_arr.sum()
    if n_pos == 0:
        raise ValueError("No positive labels found")

    tp_cumsum = np.cumsum(labels_sorted)
    fp_cumsum = np.cumsum(1 - labels_sorted)

    # Precision = TP / (TP + FP)
    precision = tp_cumsum / (tp_cumsum + fp_cumsum)

    # Recall = TP / total_positives
    recall = tp_cumsum / n_pos

    # AUPRC via trapezoidal integration
    # Prepend (0, 1) point for precision-recall curve
    recall_with_zero = np.concatenate([[0], recall])
    precision_with_one = np.concatenate([[1], precision])

    # Area under curve
    auprc = np.trapz(precision_with_one, recall_with_zero)

    return float(auprc)
