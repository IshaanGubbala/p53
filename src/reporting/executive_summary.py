"""Executive summary generator for top rescue candidates."""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import pandas as pd

from .risk_annotator import compute_risk_flags, format_risk_flags


logger = logging.getLogger(__name__)


def rank_rescues(
    df: pd.DataFrame,
    primary_key: str = "is_pareto",
    secondary_key: str = "n_rescue",
    tertiary_key: str = "ddg_gain",
) -> pd.DataFrame:
    """Rank rescue candidates by priority.

    Args:
        df: DataFrame with rescue candidates
        primary_key: Primary ranking criterion (e.g., Pareto optimality)
        secondary_key: Secondary criterion (prefer simpler rescues)
        tertiary_key: Tertiary criterion (prefer larger gain)

    Returns:
        Sorted DataFrame
    """
    df = df.copy()

    # Ensure primary key is boolean
    if primary_key in df.columns:
        df = df.sort_values(
            by=[primary_key, secondary_key, tertiary_key],
            ascending=[False, True, True],  # Pareto=True first, n_rescue ascending, ddg_gain ascending (more negative)
        )
    else:
        df = df.sort_values(
            by=[secondary_key, tertiary_key],
            ascending=[True, True],
        )

    return df.reset_index(drop=True)


def generate_executive_summary(
    candidates_df: pd.DataFrame,
    top_per_target: int = 5,
    msa_cons_map: dict[int, float] | None = None,
    dist_map: dict[tuple[int, int], float] | None = None,
    protected: set[int] | None = None,
    burial_map: dict[int, str] | None = None,
    high_cons_threshold: float = 0.8,
    near_functional_threshold: float = 10.0,
) -> pd.DataFrame:
    """Generate executive summary of top rescue candidates.

    Args:
        candidates_df: All rescue candidates
        top_per_target: Number of top rescues per target
        msa_cons_map: MSA conservation map (optional)
        dist_map: Distance map (optional)
        protected: Protected positions (optional)
        burial_map: Burial map (optional)
        high_cons_threshold: Conservation threshold
        near_functional_threshold: Distance threshold

    Returns:
        DataFrame with top rescues and risk annotations
    """
    # Compute risk flags
    annotated = compute_risk_flags(
        candidates_df,
        msa_cons_map=msa_cons_map,
        dist_map=dist_map,
        protected=protected,
        burial_map=burial_map,
        high_cons_threshold=high_cons_threshold,
        near_functional_threshold=near_functional_threshold,
    )

    # Rank within each target
    summary_rows = []

    for target in annotated["target"].unique():
        target_df = annotated[annotated["target"] == target].copy()

        # Rank rescues
        ranked = rank_rescues(target_df)

        # Take top N
        top_rescues = ranked.head(top_per_target)

        summary_rows.append(top_rescues)

    if not summary_rows:
        logger.warning("No rescue candidates found for summary")
        return pd.DataFrame()

    summary_df = pd.concat(summary_rows, ignore_index=True)

    # Add formatted risk flags
    summary_df["risk_flags"] = summary_df.apply(format_risk_flags, axis=1)

    # Select summary columns
    summary_columns = [
        "target",
        "rescue_mutations",
        "n_rescue",
        "ddg_seed",
        "ddg_gain",
        "risk",
        "risk_flags",
        "n_flags",
        "is_pareto",
    ]

    # Add multi-structure columns if present
    if "ddg_consensus" in summary_df.columns:
        summary_columns.insert(5, "ddg_consensus")
        summary_columns.insert(6, "ddg_std")

    available_columns = [col for col in summary_columns if col in summary_df.columns]

    return summary_df[available_columns]


def save_summary_csv(summary_df: pd.DataFrame, output_path: Path) -> None:
    """Save executive summary as CSV.

    Args:
        summary_df: Executive summary DataFrame
        output_path: Output CSV path
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    summary_df.to_csv(output_path, index=False)
    logger.info(f"Saved executive summary CSV: {output_path}")


def save_summary_latex(summary_df: pd.DataFrame, output_path: Path) -> None:
    """Save executive summary as LaTeX table.

    Args:
        summary_df: Executive summary DataFrame
        output_path: Output .tex path
    """
    # Create display DataFrame with formatted columns
    display_df = summary_df.copy()

    # Format numeric columns
    if "ddg_seed" in display_df.columns:
        display_df["ddg_seed"] = display_df["ddg_seed"].apply(lambda x: f"{x:+.2f}")

    if "ddg_gain" in display_df.columns:
        display_df["ddg_gain"] = display_df["ddg_gain"].apply(lambda x: f"{x:+.2f}")

    if "risk" in display_df.columns:
        display_df["risk"] = display_df["risk"].apply(lambda x: f"{x:.2f}")

    if "ddg_consensus" in display_df.columns:
        display_df["ddg_consensus"] = display_df["ddg_consensus"].apply(lambda x: f"{x:+.2f}")

    if "ddg_std" in display_df.columns:
        display_df["ddg_std"] = display_df["ddg_std"].apply(lambda x: f"{x:.2f}")

    # Rename columns for LaTeX
    column_mapping = {
        "target": "Target",
        "rescue_mutations": "Rescue Mutations",
        "n_rescue": "N",
        "ddg_seed": "$\\Delta\\Delta G_{seed}$",
        "ddg_gain": "$\\Delta\\Delta G_{gain}$",
        "ddg_consensus": "$\\Delta\\Delta G_{cons}$",
        "ddg_std": "$\\sigma$",
        "risk": "Risk",
        "risk_flags": "Flags",
        "n_flags": "$N_{flags}$",
        "is_pareto": "Pareto",
    }

    display_df = display_df.rename(columns=column_mapping)

    # Select columns for table
    table_columns = [col for col in column_mapping.values() if col in display_df.columns]
    display_df = display_df[table_columns]

    # Generate LaTeX manually (avoid jinja2 dependency)
    # Table header
    headers = " & ".join(display_df.columns) + " \\\\"
    col_format = "l" + "c" * (len(display_df.columns) - 1)

    # Table rows
    rows = []
    for _, row in display_df.iterrows():
        # Escape special LaTeX characters in string values
        escaped_values = []
        for val in row:
            if isinstance(val, str):
                val = val.replace("_", "\\_").replace("&", "\\&")
            escaped_values.append(str(val))
        rows.append(" & ".join(escaped_values) + " \\\\")

    table_content = "\n".join(rows)

    latex_table = f"""\\begin{{table}}[htbp]
\\centering
\\caption{{Executive Summary: Top Rescue Candidates}}
\\label{{tab:executive_summary}}
\\begin{{tabular}}{{{col_format}}}
\\toprule
{headers}
\\midrule
{table_content}
\\bottomrule
\\end{{tabular}}
\\end{{table}}
"""

    # Add table header
    latex_document = f"""\\documentclass{{article}}
\\usepackage{{booktabs}}
\\usepackage{{graphicx}}

\\begin{{document}}

{latex_table}

\\end{{document}}
"""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(latex_document, encoding="utf-8")
    logger.info(f"Saved executive summary LaTeX: {output_path}")
