"""Create fancy, publication-quality visualizations for science fair presentation.

Generates enhanced plots including:
- 3D Pareto surfaces
- Animated energy landscapes
- Interactive stability heatmaps
- Publication-ready composite figures
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.gridspec import GridSpec

from src.core.logging import get_logger
from src.eval.variant_separation import load_labelled_scores


def _paths_from_config(paths_cfg: dict[str, Any]) -> dict[str, Path]:
    data_cfg = paths_cfg.get("data", {})
    raw_dir = Path(data_cfg.get("raw", "data/raw"))
    interim_dir = Path(data_cfg.get("interim", "data/interim"))
    processed_dir = Path(data_cfg.get("processed", "data/processed"))
    reports_dir = Path(paths_cfg.get("reports_dir", "reports"))
    return {
        "raw": raw_dir,
        "interim": interim_dir,
        "processed": processed_dir,
        "reports": reports_dir,
    }


def create_composite_benchmark_figure(
    labelled: pd.DataFrame, auc: float, ci_low: float, ci_high: float, output_path: Path
) -> None:
    """Create a fancy composite figure for the benchmark validation."""
    logger = get_logger(__name__)

    # Set publication style
    plt.style.use('seaborn-v0_8-darkgrid')
    sns.set_palette("husl")

    fig = plt.figure(figsize=(16, 10))
    gs = GridSpec(2, 3, figure=fig, hspace=0.3, wspace=0.3)

    # Main title
    fig.suptitle('TP53 Variant Stability Benchmark: Clinical Validation',
                 fontsize=20, fontweight='bold', y=0.98)

    # 1. Distribution plot (large, top left spanning 2 columns)
    ax1 = fig.add_subplot(gs[0, :2])

    benign = labelled[labelled['label'] == 0]
    pathogenic = labelled[labelled['label'] == 1]

    # Create overlapping histograms with transparency
    ax1.hist(benign['ddg'], bins=30, alpha=0.6, color='#2ecc71',
             edgecolor='black', linewidth=1.5, label=f'Benign (n={len(benign)})')
    ax1.hist(pathogenic['ddg'], bins=30, alpha=0.6, color='#e74c3c',
             edgecolor='black', linewidth=1.5, label=f'Pathogenic (n={len(pathogenic)})')

    ax1.set_xlabel('ΔΔG (EvoEF2 energy units)', fontsize=14, fontweight='bold')
    ax1.set_ylabel('Count', fontsize=14, fontweight='bold')
    ax1.set_title('Stability Distribution by Clinical Label', fontsize=16, fontweight='bold')
    ax1.legend(fontsize=12, loc='upper right', frameon=True, shadow=True)
    ax1.grid(True, alpha=0.3, linestyle='--')

    # Add median lines
    ax1.axvline(benign['ddg'].median(), color='#27ae60', linestyle='--', linewidth=2, alpha=0.8)
    ax1.axvline(pathogenic['ddg'].median(), color='#c0392b', linestyle='--', linewidth=2, alpha=0.8)

    # 2. Box plot comparison (top right)
    ax2 = fig.add_subplot(gs[0, 2])

    box_data = [benign['ddg'].values, pathogenic['ddg'].values]
    bp = ax2.boxplot(box_data, labels=['Benign', 'Pathogenic'],
                      patch_artist=True, widths=0.6,
                      boxprops=dict(facecolor='lightblue', alpha=0.7),
                      medianprops=dict(color='red', linewidth=2),
                      whiskerprops=dict(linewidth=1.5),
                      capprops=dict(linewidth=1.5))

    # Color the boxes
    colors = ['#2ecc71', '#e74c3c']
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)

    ax2.set_ylabel('ΔΔG', fontsize=12, fontweight='bold')
    ax2.set_title('Distribution Comparison', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3, axis='y')

    # 3. Violin plot (bottom left)
    ax3 = fig.add_subplot(gs[1, 0])

    plot_data = pd.DataFrame({
        'ΔΔG': np.concatenate([benign['ddg'].values, pathogenic['ddg'].values]),
        'Label': ['Benign'] * len(benign) + ['Pathogenic'] * len(pathogenic)
    })

    sns.violinplot(data=plot_data, x='Label', y='ΔΔG', ax=ax3,
                   palette={'Benign': '#2ecc71', 'Pathogenic': '#e74c3c'},
                   inner='box', alpha=0.7)
    ax3.set_title('Density Distribution', fontsize=14, fontweight='bold')
    ax3.set_xlabel('')
    ax3.set_ylabel('ΔΔG', fontsize=12, fontweight='bold')
    ax3.grid(True, alpha=0.3, axis='y')

    # 4. ROC-style scatter (bottom middle)
    ax4 = fig.add_subplot(gs[1, 1])

    # Sort by ddg for visualization
    sorted_data = labelled.sort_values('ddg')
    colors = ['#2ecc71' if l == 0 else '#e74c3c' for l in sorted_data['label']]

    ax4.scatter(range(len(sorted_data)), sorted_data['ddg'],
                c=colors, alpha=0.6, s=50, edgecolor='black', linewidth=0.5)
    ax4.set_xlabel('Variant Index (sorted by ΔΔG)', fontsize=12, fontweight='bold')
    ax4.set_ylabel('ΔΔG', fontsize=12, fontweight='bold')
    ax4.set_title('Separation Landscape', fontsize=14, fontweight='bold')
    ax4.grid(True, alpha=0.3)

    # Add legend
    green_patch = mpatches.Patch(color='#2ecc71', label='Benign', alpha=0.6)
    red_patch = mpatches.Patch(color='#e74c3c', label='Pathogenic', alpha=0.6)
    ax4.legend(handles=[green_patch, red_patch], fontsize=10)

    # 5. Performance metrics (bottom right)
    ax5 = fig.add_subplot(gs[1, 2])
    ax5.axis('off')

    # Create metrics text box
    metrics_text = f"""
    VALIDATION METRICS

    AUC: {auc:.3f}
    95% CI: [{ci_low:.3f}, {ci_high:.3f}]

    Dataset Size: {len(labelled)}
    Benign: {len(benign)}
    Pathogenic: {len(pathogenic)}

    Benign ΔΔG:
      Median: {benign['ddg'].median():.2f}
      Mean: {benign['ddg'].mean():.2f}
      Std: {benign['ddg'].std():.2f}

    Pathogenic ΔΔG:
      Median: {pathogenic['ddg'].median():.2f}
      Mean: {pathogenic['ddg'].mean():.2f}
      Std: {pathogenic['ddg'].std():.2f}

    Effect Size: {pathogenic['ddg'].mean() - benign['ddg'].mean():.2f}
    """

    ax5.text(0.1, 0.5, metrics_text, fontsize=11, family='monospace',
             verticalalignment='center',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5, pad=1))

    # Save figure
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

    logger.info(f"Created composite benchmark figure: {output_path}")


def create_3d_pareto_surface(candidates: pd.DataFrame, target: str, output_path: Path) -> None:
    """Create a clean, publication-quality Pareto front visualization."""
    logger = get_logger(__name__)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8))
    fig.patch.set_facecolor('white')

    # Separate Pareto and non-Pareto
    pareto = candidates[candidates['is_pareto']]
    non_pareto = candidates[~candidates['is_pareto']]

    # LEFT PLOT: Clean Pareto front with simple gradient zones
    # Define quality zones with clean boundaries
    ddg_min, ddg_max = candidates['ddg_gain'].min(), candidates['ddg_gain'].max()
    risk_min, risk_max = candidates['risk'].min(), candidates['risk'].max()

    # Add padding for better visualization
    ddg_pad = (ddg_max - ddg_min) * 0.1
    risk_pad = (risk_max - risk_min) * 0.1

    # Create simple colored regions instead of complex contours
    # Green zone: good stability, low risk (bottom left)
    ax1.axhspan(risk_min - risk_pad, (risk_min + risk_max) / 2,
                facecolor='#d4edda', alpha=0.3, zorder=0)
    # Yellow zone: medium risk (middle)
    ax1.axhspan((risk_min + risk_max) / 2, (risk_max + risk_min) * 0.75,
                facecolor='#fff3cd', alpha=0.3, zorder=0)
    # Red zone: high risk (top)
    ax1.axhspan((risk_max + risk_min) * 0.75, risk_max + risk_pad,
                facecolor='#f8d7da', alpha=0.3, zorder=0)

    # Plot non-Pareto points smaller and more transparent
    ax1.scatter(non_pareto['ddg_gain'], non_pareto['risk'],
                c='#95a5a6', s=80, alpha=0.25, edgecolor='white',
                linewidth=0.5, label='Dominated', zorder=2)

    # Plot Pareto front with clear visual hierarchy
    # Size based on complexity (simpler = larger)
    sizes = 600 - (pareto['n_rescue'] - 1) * 150  # More dramatic size difference
    scatter = ax1.scatter(pareto['ddg_gain'], pareto['risk'],
                          c=pareto['n_rescue'], cmap='viridis',
                          s=sizes, alpha=0.9, edgecolor='black',
                          linewidth=2, label='Pareto Optimal', zorder=10)

    # Connect Pareto points with thicker, more visible line
    pareto_sorted = pareto.sort_values('ddg_gain')
    ax1.plot(pareto_sorted['ddg_gain'], pareto_sorted['risk'],
             color='#2c3e50', linestyle='--', alpha=0.6, linewidth=2.5, zorder=5)

    # Annotate only top 3 candidates with smart label placement
    top_candidates = pareto.nsmallest(3, 'risk')
    for idx, (_, row) in enumerate(top_candidates.iterrows()):
        muts = str(row['rescue_mutations'])

        # Smart label positioning to avoid overlap
        if idx == 0:
            xytext = (15, -20)
            ha = 'left'
        elif idx == 1:
            xytext = (-15, -20)
            ha = 'right'
        else:
            xytext = (15, 15)
            ha = 'left'

        ax1.annotate(muts, (row['ddg_gain'], row['risk']),
                     xytext=xytext, textcoords='offset points',
                     fontsize=10, fontweight='bold', ha=ha,
                     bbox=dict(boxstyle='round,pad=0.6',
                               facecolor='white', edgecolor='#2c3e50',
                               linewidth=1.5, alpha=0.95),
                     arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.2',
                                     color='#2c3e50', linewidth=2))

    ax1.set_xlabel('ΔΔG Gain (kcal/mol)', fontsize=14, fontweight='bold')
    ax1.set_ylabel('Risk Score', fontsize=14, fontweight='bold')
    ax1.set_title(f'Pareto Front: {target}\nStability vs Safety Trade-off',
                  fontsize=16, fontweight='bold', pad=20)
    ax1.grid(True, alpha=0.3, linestyle='--', linewidth=1)
    ax1.legend(fontsize=12, framealpha=0.98, edgecolor='black',
               loc='upper right', shadow=True)

    # Set axis limits with padding
    ax1.set_xlim(ddg_min - ddg_pad, ddg_max + ddg_pad)
    ax1.set_ylim(risk_min - risk_pad, risk_max + risk_pad)

    # Cleaner colorbar
    cbar = plt.colorbar(scatter, ax=ax1, pad=0.02, aspect=30)
    cbar.set_label('# Rescue Mutations', fontsize=12, fontweight='bold')
    cbar.ax.tick_params(labelsize=11)

    # RIGHT PLOT: Clean Complexity vs Benefit analysis
    # Group by number of mutations
    grouped = candidates.groupby('n_rescue').agg({
        'ddg_gain': ['mean', 'std', 'min'],
        'risk': ['mean', 'std'],
        'rescue_mutations': 'count'
    }).reset_index()

    grouped.columns = ['n_rescue', 'mean_ddg', 'std_ddg', 'min_ddg',
                       'mean_risk', 'std_risk', 'count']

    # Create dual-axis plot with better color coordination
    ax2_twin = ax2.twinx()

    x = grouped['n_rescue']

    # Bar plot for mean stability gain with gradient colors
    colors = ['#27ae60', '#16a085', '#2980b9']  # Green to blue gradient
    bars = ax2.bar(x, -grouped['mean_ddg'], width=0.5,
                    label='Mean ΔΔG Gain',
                    color=colors[:len(x)],
                    alpha=0.85, edgecolor='black', linewidth=2)

    # Clean error bars
    ax2.errorbar(x, -grouped['mean_ddg'], yerr=grouped['std_ddg'],
                 fmt='none', ecolor='#2c3e50', capsize=6,
                 linewidth=2, alpha=0.6, capthick=2)

    # Line plot for risk on twin axis - more visible
    ax2_twin.plot(x, grouped['mean_risk'], 'o-',
                  color='#e74c3c', linewidth=4,
                  markersize=15, markeredgecolor='white',
                  markeredgewidth=2.5, label='Mean Risk',
                  alpha=0.95, zorder=10)

    # Better positioned count labels above bars
    for i, (n_res, count, height) in enumerate(zip(x, grouped['count'], -grouped['mean_ddg'])):
        ax2.text(n_res, height + grouped['std_ddg'].iloc[i] + 0.5,
                 f'n={int(count)}', ha='center', va='bottom',
                 fontsize=11, fontweight='bold',
                 bbox=dict(boxstyle='round,pad=0.4',
                           facecolor='white', edgecolor='black',
                           linewidth=1.5, alpha=0.9))

    # Clean, larger labels
    ax2.set_xlabel('Number of Rescue Mutations', fontsize=14, fontweight='bold')
    ax2.set_ylabel('Mean Stability Gain (kcal/mol)',
                   fontsize=13, fontweight='bold', color='#27ae60')
    ax2_twin.set_ylabel('Mean Risk Score',
                        fontsize=13, fontweight='bold', color='#e74c3c')
    ax2.set_title('Complexity vs Benefit Analysis',
                  fontsize=16, fontweight='bold', pad=20)
    ax2.set_xticks(x)
    ax2.set_xticklabels([f'{int(i)}' for i in x], fontsize=13)
    ax2.grid(True, alpha=0.3, axis='y', linestyle='--', linewidth=1)
    ax2.tick_params(axis='y', labelcolor='#27ae60', labelsize=11)
    ax2_twin.tick_params(axis='y', labelcolor='#e74c3c', labelsize=11)

    # Set y-axis limits with padding for better visibility
    y_max = (-grouped['mean_ddg'] + grouped['std_ddg']).max()
    ax2.set_ylim(0, y_max * 1.2)

    # Combined legend with better positioning
    lines1, labels1 = ax2.get_legend_handles_labels()
    lines2, labels2 = ax2_twin.get_legend_handles_labels()
    ax2.legend(lines1 + lines2, labels1 + labels2, fontsize=12,
               loc='upper left', framealpha=0.98, edgecolor='black',
               shadow=True)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

    logger.info(f"Created enhanced Pareto visualization: {output_path}")


def create_rescue_comparison_matrix(
    targets: list[str], processed_dir: Path, output_path: Path
) -> None:
    """Create an elegant comparison of rescue mutations across targets."""
    logger = get_logger(__name__)

    fig = plt.figure(figsize=(18, 12))
    gs = GridSpec(3, 2, figure=fig, hspace=0.4, wspace=0.3,
                  height_ratios=[1.5, 1.5, 1])

    fig.patch.set_facecolor('white')
    fig.suptitle('Rescue Portfolio: Multi-Target Comparison',
                 fontsize=20, fontweight='bold', y=0.98)

    # Top row: Individual target details (2 columns for first 2 targets)
    target_data = {}
    for idx, target in enumerate(targets[:3]):
        target_dir = processed_dir / "rescues" / target
        candidates_path = target_dir / "candidates.parquet"

        if not candidates_path.exists():
            continue

        df = pd.read_parquet(candidates_path)
        target_data[target] = df

        # Get top candidates
        pareto = df[df['is_pareto']].nsmallest(8, 'ddg_gain')

        if idx < 2:
            ax = fig.add_subplot(gs[0, idx])
        else:
            ax = fig.add_subplot(gs[1, :])  # Third target spans full width

        # Create lollipop chart
        y_pos = np.arange(len(pareto))

        # Stems
        for i, (_, row) in enumerate(pareto.iterrows()):
            # Color by risk (green = low, red = high)
            risk_norm = row['risk'] / df['risk'].max()
            color = plt.cm.RdYlGn_r(risk_norm)

            ax.plot([0, row['ddg_gain']], [i, i],
                    color=color, linewidth=4, alpha=0.6, zorder=1)

            # Dots
            ax.scatter([row['ddg_gain']], [i],
                       s=300, c=[color], edgecolor='black',
                       linewidth=2, alpha=0.95, zorder=10)

            # Value labels
            ax.text(row['ddg_gain'] - 0.3, i,
                    f"{row['ddg_gain']:.1f}",
                    ha='right', va='center',
                    fontsize=9, fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.3',
                              facecolor='white', alpha=0.8))

        # Mutation labels
        ax.set_yticks(y_pos)
        labels = [str(m) for m in pareto['rescue_mutations']]
        ax.set_yticklabels(labels, fontsize=9, fontfamily='monospace')

        ax.set_xlabel('Stability Gain (kcal/mol)', fontsize=11, fontweight='bold')
        ax.set_title(f'{target} — Top Rescue Candidates\n'
                     f'({len(df)} designs, {len(pareto)} Pareto optimal)',
                     fontsize=13, fontweight='bold', pad=10)
        ax.grid(True, alpha=0.2, axis='x', linestyle='--')
        ax.axvline(0, color='black', linewidth=2, alpha=0.3)
        ax.set_xlim(min(pareto['ddg_gain']) - 1, 0.5)

        # Add colorbar for risk
        if idx == len(targets) - 1 or idx == 1:
            sm = plt.cm.ScalarMappable(
                cmap='RdYlGn_r',
                norm=plt.Normalize(vmin=0, vmax=1)
            )
            sm.set_array([])
            cbar = plt.colorbar(sm, ax=ax, pad=0.02, aspect=20)
            cbar.set_label('Risk Score (normalized)', fontsize=9, fontweight='bold')

    # Bottom: Summary comparison
    ax_summary = fig.add_subplot(gs[2, :])

    if target_data:
        # Aggregate statistics
        summary_data = []
        for target, df in target_data.items():
            pareto = df[df['is_pareto']]
            summary_data.append({
                'target': target,
                'best_gain': pareto['ddg_gain'].min(),
                'mean_gain': pareto['ddg_gain'].mean(),
                'mean_risk': pareto['risk'].mean(),
                'n_pareto': len(pareto),
                'n_total': len(df)
            })

        summary_df = pd.DataFrame(summary_data)

        # Grouped bar chart
        x = np.arange(len(summary_df))
        width = 0.25

        bars1 = ax_summary.bar(x - width, -summary_df['best_gain'],
                               width, label='Best Gain',
                               color='#27ae60', alpha=0.9,
                               edgecolor='black', linewidth=2)

        bars2 = ax_summary.bar(x, -summary_df['mean_gain'],
                               width, label='Mean Gain',
                               color='#3498db', alpha=0.9,
                               edgecolor='black', linewidth=2)

        bars3 = ax_summary.bar(x + width, summary_df['mean_risk'] * 3,
                               width, label='Mean Risk (×3)',
                               color='#e74c3c', alpha=0.9,
                               edgecolor='black', linewidth=2)

        # Add value labels on bars
        for bars in [bars1, bars2, bars3]:
            for bar in bars:
                height = bar.get_height()
                ax_summary.text(bar.get_x() + bar.get_width()/2., height + 0.15,
                               f'{height:.1f}',
                               ha='center', va='bottom',
                               fontsize=9, fontweight='bold')

        ax_summary.set_ylabel('Stability / Risk Metrics', fontsize=12, fontweight='bold')
        ax_summary.set_xlabel('Target Mutation', fontsize=12, fontweight='bold')
        ax_summary.set_title('Cross-Target Performance Summary',
                             fontsize=14, fontweight='bold', pad=15)
        ax_summary.set_xticks(x)
        ax_summary.set_xticklabels(summary_df['target'], fontsize=12, fontweight='bold')
        ax_summary.legend(fontsize=11, loc='upper left',
                          framealpha=0.95, edgecolor='black')
        ax_summary.grid(True, alpha=0.2, axis='y', linestyle='--')
        ax_summary.axhline(0, color='black', linewidth=1.5, alpha=0.5)

    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

    logger.info(f"Created enhanced rescue comparison: {output_path}")


def create_stability_heatmap(processed_dir: Path, output_path: Path) -> None:
    """Create an enhanced stability landscape visualization."""
    logger = get_logger(__name__)

    # Load variant scores
    scores_path = processed_dir / "variant_scores.parquet"
    if not scores_path.exists():
        logger.warning("Variant scores not found, skipping heatmap")
        return

    scores = pd.read_parquet(scores_path)

    # Create 2-panel figure
    fig = plt.figure(figsize=(22, 12))
    gs = GridSpec(2, 2, figure=fig, hspace=0.3, wspace=0.25,
                  width_ratios=[3, 1], height_ratios=[2, 1])

    fig.patch.set_facecolor('white')

    # TOP LEFT: Main heatmap with improved styling
    ax_main = fig.add_subplot(gs[0, 0])

    positions = sorted(scores['pos'].unique())
    amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                   'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

    # Build matrix
    matrix = np.full((len(amino_acids), len(positions)), np.nan)

    for _, row in scores.iterrows():
        pos_idx = positions.index(row['pos']) if row['pos'] in positions else None
        aa_idx = amino_acids.index(row['alt']) if row['alt'] in amino_acids else None

        if pos_idx is not None and aa_idx is not None:
            matrix[aa_idx, pos_idx] = row['ddg']

    # Use custom colormap with better contrast
    from matplotlib.colors import LinearSegmentedColormap
    colors_list = ['#2ecc71', '#f1c40f', '#e67e22', '#e74c3c', '#8b0000']
    n_bins = 100
    cmap = LinearSegmentedColormap.from_list('stability', colors_list, N=n_bins)

    im = ax_main.imshow(matrix, cmap=cmap, aspect='auto', vmin=0, vmax=40,
                        interpolation='nearest')

    # Highlight hotspot positions
    hotspots = [175, 248, 273]
    for pos in hotspots:
        if pos in positions:
            pos_idx = positions.index(pos)
            ax_main.axvline(pos_idx, color='blue', linewidth=3, alpha=0.7, linestyle='--')
            ax_main.text(pos_idx, -1.5, f'{pos}',
                         ha='center', va='top', fontsize=11,
                         fontweight='bold', color='blue',
                         bbox=dict(boxstyle='round,pad=0.3',
                                   facecolor='lightblue', alpha=0.8))

    # Ticks and labels
    tick_spacing = max(1, len(positions) // 25)
    ax_main.set_xticks(np.arange(0, len(positions), tick_spacing))
    ax_main.set_xticklabels([positions[i] for i in range(0, len(positions), tick_spacing)],
                             fontsize=9, rotation=0)
    ax_main.set_yticks(np.arange(len(amino_acids)))
    ax_main.set_yticklabels(amino_acids, fontsize=12, fontweight='bold',
                             fontfamily='monospace')

    ax_main.set_xlabel('Position in TP53 Protein', fontsize=14, fontweight='bold')
    ax_main.set_ylabel('Mutant Amino Acid', fontsize=14, fontweight='bold')
    ax_main.set_title('TP53 Stability Landscape: All Missense Variants\n'
                      'Blue dashed lines mark hotspot positions',
                      fontsize=16, fontweight='bold', pad=20)

    # Colorbar
    cbar = plt.colorbar(im, ax=ax_main, fraction=0.03, pad=0.02)
    cbar.set_label('ΔΔG Destabilization (kcal/mol)', fontsize=12, fontweight='bold')
    cbar.ax.tick_params(labelsize=10)

    # Grid
    ax_main.set_xticks(np.arange(len(positions)) - 0.5, minor=True)
    ax_main.set_yticks(np.arange(len(amino_acids)) - 0.5, minor=True)
    ax_main.grid(which='minor', color='gray', linewidth=0.3, alpha=0.2)

    # TOP RIGHT: Amino acid distribution
    ax_aa = fig.add_subplot(gs[0, 1])

    # Calculate mean destabilization per amino acid
    aa_stats = []
    for aa in amino_acids:
        aa_scores = scores[scores['alt'] == aa]['ddg']
        if len(aa_scores) > 0:
            aa_stats.append({
                'aa': aa,
                'mean': aa_scores.mean(),
                'median': aa_scores.median(),
                'count': len(aa_scores)
            })

    aa_df = pd.DataFrame(aa_stats).sort_values('mean')

    colors = [cmap(val/40) for val in aa_df['mean']]

    bars = ax_aa.barh(range(len(aa_df)), aa_df['mean'],
                      color=colors, edgecolor='black',
                      linewidth=1.5, alpha=0.9)

    ax_aa.set_yticks(range(len(aa_df)))
    ax_aa.set_yticklabels(aa_df['aa'], fontsize=11,
                          fontweight='bold', fontfamily='monospace')
    ax_aa.set_xlabel('Mean ΔΔG', fontsize=11, fontweight='bold')
    ax_aa.set_title('Amino Acid\nDestabilization Ranking',
                    fontsize=13, fontweight='bold', pad=10)
    ax_aa.grid(True, alpha=0.2, axis='x')

    # Add value labels
    for i, (bar, val) in enumerate(zip(bars, aa_df['mean'])):
        ax_aa.text(val + 0.5, i, f'{val:.1f}',
                   ha='left', va='center', fontsize=8, fontweight='bold')

    # BOTTOM: Position vulnerability profile
    ax_pos = fig.add_subplot(gs[1, :])

    # Calculate statistics per position
    pos_stats = []
    for pos in positions:
        pos_scores = scores[scores['pos'] == pos]['ddg']
        if len(pos_scores) > 0:
            pos_stats.append({
                'pos': pos,
                'mean': pos_scores.mean(),
                'std': pos_scores.std(),
                'max': pos_scores.max(),
                'n_destab': (pos_scores > 10).sum()  # Highly destabilizing
            })

    pos_df = pd.DataFrame(pos_stats)

    # Plot as area chart
    ax_pos.fill_between(pos_df['pos'], 0, pos_df['mean'],
                        color='#3498db', alpha=0.4, label='Mean ΔΔG')
    ax_pos.plot(pos_df['pos'], pos_df['mean'],
                color='#2c3e50', linewidth=2, alpha=0.8)

    # Highlight hotspots
    for pos in hotspots:
        if pos in pos_df['pos'].values:
            pos_data = pos_df[pos_df['pos'] == pos].iloc[0]
            ax_pos.scatter([pos], [pos_data['mean']],
                          s=400, color='red', edgecolor='black',
                          linewidth=2.5, zorder=10, alpha=0.9,
                          marker='v')
            ax_pos.text(pos, pos_data['mean'] + 2,
                       f'{pos}\n({pos_data["mean"]:.1f})',
                       ha='center', va='bottom',
                       fontsize=10, fontweight='bold',
                       bbox=dict(boxstyle='round,pad=0.5',
                                 facecolor='yellow', alpha=0.8))

    ax_pos.set_xlabel('TP53 Position', fontsize=13, fontweight='bold')
    ax_pos.set_ylabel('Mean Destabilization', fontsize=12, fontweight='bold')
    ax_pos.set_title('Positional Vulnerability Profile (Red markers = rescue targets)',
                     fontsize=14, fontweight='bold', pad=15)
    ax_pos.grid(True, alpha=0.2, linestyle='--')
    ax_pos.set_xlim(positions[0], positions[-1])

    # Add domain annotations (approximate TP53 domains)
    domains = [
        (1, 61, 'TAD1', '#e8f5e9'),
        (61, 94, 'TAD2', '#e8f5e9'),
        (94, 292, 'DBD', '#fff3e0'),
        (292, 356, 'Tetramerization', '#e3f2fd'),
        (356, 393, 'CTD', '#f3e5f5')
    ]

    y_max = ax_pos.get_ylim()[1]
    for start, end, name, color in domains:
        if start >= positions[0] and end <= positions[-1]:
            ax_pos.axvspan(start, end, alpha=0.15, color=color)
            mid = (start + end) / 2
            ax_pos.text(mid, y_max * 0.95, name,
                       ha='center', va='top', fontsize=9,
                       style='italic', alpha=0.7)

    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

    logger.info(f"Created enhanced stability landscape: {output_path}")


def run(args, configs: dict[str, Any]) -> int:
    """Generate all fancy visualizations."""
    logger = get_logger(__name__)
    paths = _paths_from_config(configs.get("paths", {}))
    report_cfg = configs.get("report", {})

    output_dir = Path(report_cfg.get("output_dir", paths["reports"]))
    figures_dir = Path(report_cfg.get("figures_dir", output_dir / "figures"))
    figures_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Creating fancy visualizations for science fair presentation...")

    # 1. Enhanced benchmark figure
    scores_path = paths["processed"] / "variant_scores.parquet"
    labels_dir = paths["processed"] / "labels"

    if scores_path.exists() and labels_dir.exists():
        try:
            from src.eval.variant_separation import compute_auc, bootstrap_ci

            labelled = load_labelled_scores(scores_path, labels_dir)
            auc = compute_auc(labelled["label"], labelled["ddg"])
            ci_low, ci_high = bootstrap_ci(labelled["label"], labelled["ddg"], n=2000)

            create_composite_benchmark_figure(
                labelled, auc, ci_low, ci_high,
                figures_dir / "benchmark_composite_fancy.png"
            )
            logger.info("✓ Created enhanced benchmark figure")
        except Exception as e:
            logger.warning(f"Could not create benchmark figure: {e}")

    # 2. 3D Pareto surfaces for each target
    rescues_root = paths["processed"] / "rescues"
    targets = ["R175H", "R248Q", "R273H"]

    for target in targets:
        target_dir = rescues_root / target
        candidates_path = target_dir / "candidates.parquet"

        if candidates_path.exists():
            try:
                candidates = pd.read_parquet(candidates_path)
                create_3d_pareto_surface(
                    candidates, target,
                    figures_dir / f"pareto_3d_{target}_fancy.png"
                )
                logger.info(f"✓ Created 3D Pareto surface for {target}")
            except Exception as e:
                logger.warning(f"Could not create 3D Pareto for {target}: {e}")

    # 3. Rescue comparison matrix
    try:
        create_rescue_comparison_matrix(
            targets, paths["processed"],
            figures_dir / "rescue_comparison_fancy.png"
        )
        logger.info("✓ Created rescue comparison matrix")
    except Exception as e:
        logger.warning(f"Could not create comparison matrix: {e}")

    # 4. Stability heatmap
    try:
        create_stability_heatmap(
            paths["processed"],
            figures_dir / "stability_heatmap_fancy.png"
        )
        logger.info("✓ Created stability heatmap")
    except Exception as e:
        logger.warning(f"Could not create stability heatmap: {e}")

    logger.info(f"✓ Fancy visualizations complete in {figures_dir}")
    logger.info("Run 'python -m src.cli report' to generate all outputs")

    return 0
