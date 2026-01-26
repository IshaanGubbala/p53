"""Train/test splitting for ClinVar variants with stratification."""
from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Optional
import numpy as np
import pandas as pd


logger = logging.getLogger(__name__)


def stratified_train_test_split(
    df: pd.DataFrame,
    stratify_column: str,
    test_ratio: float = 0.2,
    random_seed: int = 42,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Perform stratified train/test split.

    Args:
        df: DataFrame with variants
        stratify_column: Column to stratify by (e.g., 'clinical_significance')
        test_ratio: Fraction of data for test set (0-1)
        random_seed: Random seed for reproducibility

    Returns:
        Tuple of (train_df, test_df)
    """
    np.random.seed(random_seed)

    train_indices = []
    test_indices = []

    # Split within each stratum
    for stratum_value in df[stratify_column].unique():
        stratum_df = df[df[stratify_column] == stratum_value]
        n_total = len(stratum_df)
        n_test = max(1, int(n_total * test_ratio))  # At least 1 per stratum

        # Shuffle indices
        indices = stratum_df.index.tolist()
        np.random.shuffle(indices)

        # Split
        test_indices.extend(indices[:n_test])
        train_indices.extend(indices[n_test:])

    train_df = df.loc[train_indices]
    test_df = df.loc[test_indices]

    logger.info(f"Stratified split by '{stratify_column}':")
    logger.info(f"  Train: {len(train_df)} ({100 * len(train_df) / len(df):.1f}%)")
    logger.info(f"  Test:  {len(test_df)} ({100 * len(test_df) / len(df):.1f}%)")

    # Log stratification balance
    for stratum_value in df[stratify_column].unique():
        train_count = (train_df[stratify_column] == stratum_value).sum()
        test_count = (test_df[stratify_column] == stratum_value).sum()
        total_count = (df[stratify_column] == stratum_value).sum()

        logger.info(
            f"  {stratum_value}: train={train_count}/{total_count} "
            f"({100 * train_count / total_count:.1f}%), "
            f"test={test_count}/{total_count} ({100 * test_count / total_count:.1f}%)"
        )

    return train_df, test_df


def create_split_indices(
    benign_path: Path,
    pathogenic_path: Path,
    split_ratio: float = 0.8,
    split_seed: int = 42,
    stratify_by: Optional[str] = None,
) -> dict[str, any]:
    """Create train/test split indices for ClinVar labels.

    Args:
        benign_path: Path to benign.parquet
        pathogenic_path: Path to pathogenic.parquet
        split_ratio: Train fraction (0-1)
        split_seed: Random seed
        stratify_by: Column to stratify by (optional)

    Returns:
        Dictionary with split indices and metadata
    """
    logger.info("Creating train/test split...")

    # Load labels
    benign = pd.read_parquet(benign_path)
    pathogenic = pd.read_parquet(pathogenic_path)

    # Add label column
    benign["label"] = 0
    pathogenic["label"] = 1

    # Combine
    combined = pd.concat([benign, pathogenic], ignore_index=True)

    logger.info(f"Total variants: {len(combined)}")
    logger.info(f"  Benign: {len(benign)}")
    logger.info(f"  Pathogenic: {len(pathogenic)}")

    # Perform split
    if stratify_by and stratify_by in combined.columns:
        train_df, test_df = stratified_train_test_split(
            combined,
            stratify_by,
            test_ratio=1.0 - split_ratio,
            random_seed=split_seed,
        )
    else:
        # Simple random split
        np.random.seed(split_seed)
        indices = np.arange(len(combined))
        np.random.shuffle(indices)

        n_train = int(len(combined) * split_ratio)
        train_indices = indices[:n_train]
        test_indices = indices[n_train:]

        train_df = combined.iloc[train_indices]
        test_df = combined.iloc[test_indices]

        logger.info(f"Random split (no stratification):")
        logger.info(f"  Train: {len(train_df)} ({100 * len(train_df) / len(combined):.1f}%)")
        logger.info(f"  Test:  {len(test_df)} ({100 * len(test_df) / len(combined):.1f}%)")

    # Create split info
    split_info = {
        "split_ratio": split_ratio,
        "split_seed": split_seed,
        "stratify_by": stratify_by,
        "n_total": len(combined),
        "n_train": len(train_df),
        "n_test": len(test_df),
        "train_indices": train_df.index.tolist(),
        "test_indices": test_df.index.tolist(),
        "train_benign": (train_df["label"] == 0).sum(),
        "train_pathogenic": (train_df["label"] == 1).sum(),
        "test_benign": (test_df["label"] == 0).sum(),
        "test_pathogenic": (test_df["label"] == 1).sum(),
    }

    return split_info


def save_split(split_info: dict[str, any], output_path: Path) -> None:
    """Save split indices to JSON.

    Args:
        split_info: Split information dictionary
        output_path: Output JSON path
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Convert numpy types to Python native types for JSON serialization
    def convert_to_native(obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, dict):
            return {k: convert_to_native(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [convert_to_native(item) for item in obj]
        return obj

    split_info_native = convert_to_native(split_info)
    output_path.write_text(json.dumps(split_info_native, indent=2), encoding="utf-8")

    logger.info(f"Saved split to {output_path}")


def load_split(split_path: Path) -> dict[str, any]:
    """Load split indices from JSON.

    Args:
        split_path: Path to split JSON file

    Returns:
        Split information dictionary
    """
    if not split_path.exists():
        raise FileNotFoundError(f"Split file not found: {split_path}")

    split_info = json.loads(split_path.read_text(encoding="utf-8"))

    logger.info(f"Loaded split from {split_path}")
    logger.info(f"  Train: {split_info['n_train']} variants")
    logger.info(f"  Test:  {split_info['n_test']} variants")

    return split_info


def apply_split(
    df: pd.DataFrame,
    split_info: dict[str, any],
    subset: str = "train",
) -> pd.DataFrame:
    """Apply split to dataframe.

    Args:
        df: DataFrame to split
        split_info: Split information from load_split()
        subset: "train" or "test"

    Returns:
        Subset dataframe
    """
    if subset == "train":
        indices = split_info["train_indices"]
    elif subset == "test":
        indices = split_info["test_indices"]
    else:
        raise ValueError(f"Invalid subset: {subset} (must be 'train' or 'test')")

    # Filter to indices that exist in df
    valid_indices = [i for i in indices if i in df.index]

    if len(valid_indices) < len(indices):
        logger.warning(
            f"Only {len(valid_indices)}/{len(indices)} {subset} indices found in dataframe"
        )

    return df.loc[valid_indices]
