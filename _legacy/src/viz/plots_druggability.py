"""
Visualization functions for druggability analysis.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from src.core.logging import get_logger

logger = get_logger(__name__)


def plot_druggability_distribution(
    druggability_df: pd.DataFrame,
    output_path: Path,
    target: str | None = None,
    figsize: tuple[float, float] = (10, 6),
) -> None:
    """
    Plot distribution of druggability scores.

    Args:
        druggability_df: DataFrame with druggability_score column
        output_path: Path to save figure
        target: Target name for title (optional)
        figsize: Figure size
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

    # Left: Histogram
    scores = druggability_df["druggability_score"].dropna()

    ax1.hist(scores, bins=20, color='skyblue', edgecolor='black', alpha=0.7)
    ax1.axvline(0.6, color='green', linestyle='--', linewidth=2, label='High (≥0.6)')
    ax1.axvline(0.4, color='orange', linestyle='--', linewidth=2, label='Medium (≥0.4)')
    ax1.set_xlabel('Druggability Score', fontsize=12)
    ax1.set_ylabel('Count', fontsize=12)
    ax1.set_title(f'Druggability Score Distribution{f" - {target}" if target else ""}', fontsize=14)
    ax1.legend()
    ax1.grid(axis='y', alpha=0.3)

    # Right: Scatter plot (ΔΔG vs Druggability)
    if 'ddg_gain' in druggability_df.columns:
        x = druggability_df['ddg_gain']
        y = druggability_df['druggability_score']

        # Color by classification
        colors = []
        for score in y:
            if pd.isna(score):
                colors.append('gray')
            elif score >= 0.6:
                colors.append('green')
            elif score >= 0.4:
                colors.append('orange')
            else:
                colors.append('red')

        ax2.scatter(x, y, c=colors, s=100, alpha=0.6, edgecolors='black')
        ax2.axhline(0.6, color='green', linestyle='--', linewidth=1, alpha=0.5)
        ax2.axhline(0.4, color='orange', linestyle='--', linewidth=1, alpha=0.5)
        ax2.set_xlabel('ΔΔG Gain (kcal/mol)', fontsize=12)
        ax2.set_ylabel('Druggability Score', fontsize=12)
        ax2.set_title('Stability vs Druggability', fontsize=14)
        ax2.grid(alpha=0.3)

        # Add legend
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='green', label='High (≥0.6)'),
            Patch(facecolor='orange', label='Medium (0.4-0.6)'),
            Patch(facecolor='red', label='Low (<0.4)')
        ]
        ax2.legend(handles=legend_elements, loc='best')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"Druggability distribution plot saved to {output_path}")


def plot_druggability_comparison(
    all_targets_df: pd.DataFrame,
    output_path: Path,
    figsize: tuple[float, float] = (12, 6),
) -> None:
    """
    Compare druggability across targets.

    Args:
        all_targets_df: DataFrame with target and druggability_score columns
        output_path: Path to save figure
        figsize: Figure size
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

    # Left: Box plot by target
    targets = sorted(all_targets_df['target'].unique())

    data_by_target = [
        all_targets_df[all_targets_df['target'] == t]['druggability_score'].dropna()
        for t in targets
    ]

    bp = ax1.boxplot(data_by_target, labels=targets, patch_artist=True)

    for patch in bp['boxes']:
        patch.set_facecolor('lightblue')
        patch.set_alpha(0.7)

    ax1.axhline(0.6, color='green', linestyle='--', linewidth=1, alpha=0.5, label='High threshold')
    ax1.axhline(0.4, color='orange', linestyle='--', linewidth=1, alpha=0.5, label='Medium threshold')
    ax1.set_ylabel('Druggability Score', fontsize=12)
    ax1.set_xlabel('Target', fontsize=12)
    ax1.set_title('Druggability by Target', fontsize=14)
    ax1.legend()
    ax1.grid(axis='y', alpha=0.3)

    # Right: Counts by classification
    classification_counts = []
    for target in targets:
        target_df = all_targets_df[all_targets_df['target'] == target]
        scores = target_df['druggability_score'].dropna()

        high = (scores >= 0.6).sum()
        medium = ((scores >= 0.4) & (scores < 0.6)).sum()
        low = (scores < 0.4).sum()

        classification_counts.append({
            'Target': target,
            'High (≥0.6)': high,
            'Medium (0.4-0.6)': medium,
            'Low (<0.4)': low,
        })

    counts_df = pd.DataFrame(classification_counts)

    x = np.arange(len(targets))
    width = 0.25

    ax2.bar(x - width, counts_df['High (≥0.6)'], width, label='High (≥0.6)', color='green', alpha=0.7)
    ax2.bar(x, counts_df['Medium (0.4-0.6)'], width, label='Medium (0.4-0.6)', color='orange', alpha=0.7)
    ax2.bar(x + width, counts_df['Low (<0.4)'], width, label='Low (<0.4)', color='red', alpha=0.7)

    ax2.set_xlabel('Target', fontsize=12)
    ax2.set_ylabel('Count', fontsize=12)
    ax2.set_title('Classification Counts by Target', fontsize=14)
    ax2.set_xticks(x)
    ax2.set_xticklabels(targets)
    ax2.legend()
    ax2.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"Druggability comparison plot saved to {output_path}")


def plot_top_druggable_rescues(
    druggability_df: pd.DataFrame,
    output_path: Path,
    target: str | None = None,
    top_n: int = 10,
    figsize: tuple[float, float] = (10, 8),
) -> None:
    """
    Plot top druggable rescues with multiple metrics.

    Args:
        druggability_df: DataFrame with rescue data
        output_path: Path to save figure
        target: Target name for title
        top_n: Number of top rescues to show
        figsize: Figure size
    """
    # Sort by druggability score
    top_df = druggability_df.nlargest(top_n, 'druggability_score').copy()

    # Reverse order for bottom-up plotting
    top_df = top_df.iloc[::-1]

    fig, ax = plt.subplots(figsize=figsize)

    y_pos = np.arange(len(top_df))

    # Main bars (druggability score)
    bars = ax.barh(y_pos, top_df['druggability_score'], alpha=0.7, color='skyblue', edgecolor='black')

    # Color bars by classification
    for i, (idx, row) in enumerate(top_df.iterrows()):
        score = row['druggability_score']
        if score >= 0.6:
            bars[i].set_color('green')
        elif score >= 0.4:
            bars[i].set_color('orange')
        else:
            bars[i].set_color('red')

    # Add vertical lines for thresholds
    ax.axvline(0.6, color='green', linestyle='--', linewidth=2, alpha=0.5, label='High (≥0.6)')
    ax.axvline(0.4, color='orange', linestyle='--', linewidth=2, alpha=0.5, label='Medium (≥0.4)')

    # Add score labels
    for i, (idx, row) in enumerate(top_df.iterrows()):
        score = row['druggability_score']
        ax.text(score + 0.01, i, f'{score:.3f}', va='center', fontsize=10)

    # Set labels
    rescue_labels = top_df['rescue_mutations'].tolist()
    ax.set_yticks(y_pos)
    ax.set_yticklabels(rescue_labels, fontsize=10)
    ax.set_xlabel('Druggability Score', fontsize=12)
    ax.set_ylabel('Rescue Mutation', fontsize=12)
    ax.set_title(f'Top {top_n} Druggable Rescues{f" - {target}" if target else ""}', fontsize=14)
    ax.legend(loc='lower right')
    ax.grid(axis='x', alpha=0.3)
    ax.set_xlim(0, 1.0)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"Top druggable rescues plot saved to {output_path}")
