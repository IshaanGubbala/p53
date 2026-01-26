"""Risk flag annotation for rescue mutation candidates."""
from __future__ import annotations

import re
from typing import Any

import pandas as pd


def _parse_mutation(mut: str) -> tuple[int, str, str]:
    """Parse mutation string into components."""
    match = re.match(r"^([A-Za-z\*])(\d+)([A-Za-z\*])$", mut.strip())
    if not match:
        raise ValueError(f"Invalid mutation: {mut}")
    ref, pos_text, alt = match.groups()
    return int(pos_text), ref.upper(), alt.upper()


def annotate_high_conservation(
    df: pd.DataFrame,
    msa_cons_map: dict[int, float],
    threshold: float = 0.8,
) -> pd.DataFrame:
    """Annotate rescues with high MSA conservation flag.

    Args:
        df: DataFrame with rescue_mutations column
        msa_cons_map: Position → conservation score
        threshold: Conservation threshold for flagging

    Returns:
        DataFrame with added 'high_conservation' column
    """
    df = df.copy()

    def check_conservation(rescue_str: str) -> bool:
        if not rescue_str or pd.isna(rescue_str):
            return False
        mutations = rescue_str.split(",")
        for mut in mutations:
            try:
                pos, _, _ = _parse_mutation(mut)
                cons = msa_cons_map.get(pos, 0.0)
                if cons > threshold:
                    return True
            except ValueError:
                continue
        return False

    df["high_conservation"] = df["rescue_mutations"].apply(check_conservation)
    return df


def annotate_near_functional(
    df: pd.DataFrame,
    dist_map: dict[tuple[int, int], float],
    protected: set[int],
    threshold: float = 10.0,
) -> pd.DataFrame:
    """Annotate rescues near functional sites.

    Args:
        df: DataFrame with rescue_mutations column
        dist_map: Distance map (pos_i, pos_j) → distance (Angstroms)
        protected: Set of protected functional positions
        threshold: Distance threshold in Angstroms

    Returns:
        DataFrame with added 'near_functional' column
    """
    df = df.copy()

    def check_proximity(rescue_str: str) -> bool:
        if not rescue_str or pd.isna(rescue_str):
            return False
        mutations = rescue_str.split(",")
        for mut in mutations:
            try:
                pos, _, _ = _parse_mutation(mut)
                # Check distance to all protected positions
                for prot_pos in protected:
                    dist = dist_map.get((pos, prot_pos))
                    if dist is None:
                        dist = dist_map.get((prot_pos, pos))
                    if dist is not None and dist < threshold:
                        return True
            except ValueError:
                continue
        return False

    df["near_functional"] = df["rescue_mutations"].apply(check_proximity)
    return df


def annotate_exposed_surface(
    df: pd.DataFrame,
    burial_map: dict[int, str],
) -> pd.DataFrame:
    """Annotate rescues with exposed surface mutations.

    Args:
        df: DataFrame with rescue_mutations column
        burial_map: Position → burial category (buried/partial/exposed)

    Returns:
        DataFrame with added 'exposed' column
    """
    df = df.copy()

    def check_exposure(rescue_str: str) -> bool:
        if not rescue_str or pd.isna(rescue_str):
            return False
        mutations = rescue_str.split(",")
        for mut in mutations:
            try:
                pos, _, _ = _parse_mutation(mut)
                burial = burial_map.get(pos, "unknown")
                if burial == "exposed":
                    return True
            except ValueError:
                continue
        return False

    df["exposed"] = df["rescue_mutations"].apply(check_exposure)
    return df


def compute_risk_flags(
    df: pd.DataFrame,
    msa_cons_map: dict[int, float] | None = None,
    dist_map: dict[tuple[int, int], float] | None = None,
    protected: set[int] | None = None,
    burial_map: dict[int, str] | None = None,
    high_cons_threshold: float = 0.8,
    near_functional_threshold: float = 10.0,
) -> pd.DataFrame:
    """Compute all risk flags for rescue candidates.

    Args:
        df: DataFrame with rescue_mutations column
        msa_cons_map: MSA conservation scores (optional)
        dist_map: Distance map (optional)
        protected: Protected positions (optional)
        burial_map: Burial categories (optional)
        high_cons_threshold: Conservation threshold
        near_functional_threshold: Distance threshold (Angstroms)

    Returns:
        DataFrame with risk flag columns added
    """
    result = df.copy()

    # Annotate conservation
    if msa_cons_map:
        result = annotate_high_conservation(
            result, msa_cons_map, threshold=high_cons_threshold
        )
    else:
        result["high_conservation"] = False

    # Annotate functional proximity
    if dist_map and protected:
        result = annotate_near_functional(
            result, dist_map, protected, threshold=near_functional_threshold
        )
    else:
        result["near_functional"] = False

    # Annotate surface exposure
    if burial_map:
        result = annotate_exposed_surface(result, burial_map)
    else:
        result["exposed"] = False

    # Count total flags
    result["n_flags"] = (
        result["high_conservation"].astype(int)
        + result["near_functional"].astype(int)
        + result["exposed"].astype(int)
    )

    return result


def format_risk_flags(row: pd.Series) -> str:
    """Format risk flags as human-readable string.

    Args:
        row: DataFrame row with flag columns

    Returns:
        Comma-separated flag descriptions
    """
    flags = []

    if row.get("high_conservation", False):
        flags.append("High conservation")

    if row.get("near_functional", False):
        flags.append("Near functional site")

    if row.get("exposed", False):
        flags.append("Surface exposed")

    if not flags:
        return "None"

    return ", ".join(flags)
