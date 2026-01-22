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
