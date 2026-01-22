from __future__ import annotations

from pathlib import Path

import pandas as pd

from src.variants.normalize_variants import _review_status_score


BENIGN_TERMS = {
    "benign",
    "likely benign",
    "benign/likely benign",
}
PATHOGENIC_TERMS = {
    "pathogenic",
    "likely pathogenic",
    "pathogenic/likely pathogenic",
}


def _match_terms(series: pd.Series, terms: set[str]) -> pd.Series:
    cleaned = series.fillna("").str.lower()
    return cleaned.isin(terms)


def make_benign_set(df: pd.DataFrame) -> pd.DataFrame:
    if "clinical_significance" not in df.columns:
        return df.head(0).copy()
    mask = _match_terms(df["clinical_significance"], BENIGN_TERMS)
    return df.loc[mask].copy()


def make_pathogenic_set(df: pd.DataFrame) -> pd.DataFrame:
    if "clinical_significance" not in df.columns:
        return df.head(0).copy()
    mask = _match_terms(df["clinical_significance"], PATHOGENIC_TERMS)
    return df.loc[mask].copy()


def quality_filter(df: pd.DataFrame, min_review_stars: int = 1) -> pd.DataFrame:
    if "review_status" not in df.columns:
        return df.copy()

    scores = df["review_status"].apply(_review_status_score)
    return df.loc[scores >= min_review_stars].copy()


def export_label_sets(benign_df: pd.DataFrame, path_df: pd.DataFrame, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    benign_df.to_parquet(out_dir / "benign.parquet", index=False)
    path_df.to_parquet(out_dir / "pathogenic.parquet", index=False)
