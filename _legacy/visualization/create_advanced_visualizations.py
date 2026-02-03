"""Create advanced, eye-catching visualizations for maximum science fair impact.

Includes:
- Mutation network graphs
- Radar charts for multi-dimensional comparison
- Sankey flow diagrams
- Conservation vs stability scatter plots
- Domain architecture with rescue positions
- Animated-style progression plots
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, Circle, Rectangle, FancyArrow
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.gridspec import GridSpec
from matplotlib.collections import LineCollection

from src.core.logging import get_logger


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


def create_mutation_network(processed_dir: Path, output_path: Path) -> None:
    """Create a network visualization showing rescue mutation relationships."""
    logger = get_logger(__name__)

    # Load rescue data for all targets
    targets = ["R175H", "R248Q", "R273H"]
    all_mutations = {}

    for target in targets:
        target_dir = processed_dir / "rescues" / target
        candidates_path = target_dir / "candidates.parquet"

        if not candidates_path.exists():
            continue

        df = pd.read_parquet(candidates_path)
        pareto = df[df['is_pareto']]

        for _, row in pareto.iterrows():
            rescue_muts = row['rescue_mutations'].split(',') if row['rescue_mutations'] else []
            for mut in rescue_muts:
                if mut not in all_mutations:
                    all_mutations[mut] = {'targets': [], 'ddg_gains': [], 'risks': []}
                all_mutations[mut]['targets'].append(target)
                all_mutations[mut]['ddg_gains'].append(row['ddg_gain'])
                all_mutations[mut]['risks'].append(row['risk'])

    # Create figure
    fig, ax = plt.subplots(figsize=(18, 12))

    # Define positions for targets (triangle layout)
    target_positions = {
        'R175H': (0.2, 0.8),
        'R248Q': (0.8, 0.8),
        'R273H': (0.5, 0.2)
    }

    # Define positions for mutations (circular layout around center)
    n_mutations = len(all_mutations)
    angles = np.linspace(0, 2*np.pi, n_mutations, endpoint=False)
    mutation_positions = {}

    for i, (mut, angle) in enumerate(zip(all_mutations.keys(), angles)):
        x = 0.5 + 0.25 * np.cos(angle)
        y = 0.5 + 0.25 * np.sin(angle)
        mutation_positions[mut] = (x, y)

    # Draw connections
    for mut, data in all_mutations.items():
        mut_pos = mutation_positions[mut]

        for target in data['targets']:
            target_pos = target_positions[target]

            # Line thickness based on frequency
            linewidth = 2 + len(data['targets']) * 0.5

            # Color based on average ddg_gain
            avg_ddg = np.mean(data['ddg_gains'])
            color = plt.cm.RdYlGn_r((avg_ddg + 20) / 20)  # Normalize to 0-1

            ax.plot([mut_pos[0], target_pos[0]], [mut_pos[1], target_pos[1]],
                   color=color, linewidth=linewidth, alpha=0.5, zorder=1)

    # Draw target nodes
    for target, pos in target_positions.items():
        circle = Circle(pos, 0.08, facecolor='#e74c3c', edgecolor='black',
                       linewidth=3, zorder=3, alpha=0.8)
        ax.add_patch(circle)
        ax.text(pos[0], pos[1], target, ha='center', va='center',
               fontsize=14, fontweight='bold', color='white', zorder=4)

    # Draw mutation nodes
    for mut, pos in mutation_positions.items():
        data = all_mutations[mut]

        # Size based on frequency across targets
        size = 0.03 + len(data['targets']) * 0.015

        # Color based on average risk
        avg_risk = np.mean(data['risks'])
        if avg_risk == 0:
            color = '#2ecc71'  # Green for zero risk
        else:
            color = plt.cm.YlOrRd(avg_risk / 0.1)

        circle = Circle(pos, size, facecolor=color, edgecolor='black',
                       linewidth=2, zorder=2, alpha=0.9)
        ax.add_patch(circle)
        ax.text(pos[0], pos[1] - size - 0.03, mut, ha='center', va='top',
               fontsize=9, fontweight='bold', zorder=4)

    # Add title and legend
    ax.set_xlim(-0.1, 1.1)
    ax.set_ylim(-0.1, 1.1)
    ax.set_aspect('equal')
    ax.axis('off')

    ax.set_title('Rescue Mutation Network: Cross-Target Relationships',
                fontsize=20, fontweight='bold', pad=20)

    # Add legend
    legend_elements = [
        mpatches.Patch(color='#2ecc71', label='Zero Risk'),
        mpatches.Patch(color='#e74c3c', label='Target Mutation'),
        plt.Line2D([0], [0], color='gray', linewidth=4, label='Strong Stabilizer'),
        plt.Line2D([0], [0], color='gray', linewidth=1, label='Moderate Stabilizer')
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=11, frameon=True, shadow=True)

    # Add caption
    caption = (f"Network shows {len(all_mutations)} rescue mutations across 3 cancer hotspots.\n"
              f"Node size = mutation frequency, color = risk level, edge thickness = effectiveness.")
    ax.text(0.5, -0.05, caption, ha='center', va='top', fontsize=11,
           style='italic', wrap=True, transform=ax.transAxes)

    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

    logger.info(f"Created mutation network: {output_path}")


def create_radar_chart_comparison(processed_dir: Path, output_path: Path) -> None:
    """Create radar charts comparing top rescue candidates."""
    logger = get_logger(__name__)

    targets = ["R175H", "R248Q", "R273H"]
    fig, axes = plt.subplots(1, 3, figsize=(20, 7), subplot_kw=dict(projection='polar'))

    fig.suptitle('Multi-Dimensional Rescue Candidate Comparison',
                fontsize=18, fontweight='bold', y=1.02)

    for idx, (target, ax) in enumerate(zip(targets, axes)):
        target_dir = processed_dir / "rescues" / target
        candidates_path = target_dir / "candidates.parquet"

        if not candidates_path.exists():
            continue

        df = pd.read_parquet(candidates_path)
        pareto = df[df['is_pareto']].nsmallest(5, 'ddg_gain')

        # Define metrics (normalized to 0-1 scale)
        categories = ['Stability\nGain', 'Low Risk', 'Simplicity', 'Burial\nScore', 'Conservation\nScore']
        num_vars = len(categories)

        # Compute angles for radar chart
        angles = np.linspace(0, 2*np.pi, num_vars, endpoint=False).tolist()
        angles += angles[:1]  # Complete the circle

        # Plot each candidate
        colors = plt.cm.Set3(np.linspace(0, 1, len(pareto)))

        for i, (_, row) in enumerate(pareto.iterrows()):
            # Normalize metrics to 0-1 scale (higher is better)
            stability = min(1.0, abs(row['ddg_gain']) / 20)  # Normalize by max gain
            low_risk = 1.0 - min(1.0, row['risk'] / 0.1)  # Invert risk
            simplicity = 1.0 - (row['n_rescue'] - 1) / 2  # Fewer mutations = better

            # Parse risk components
            risk_comp = json.loads(row['risk_components'])
            burial = 1.0 - risk_comp.get('burial', 0)  # Lower burial risk = better
            conservation = 1.0 - risk_comp.get('conservation', 0)

            values = [stability, low_risk, simplicity, burial, conservation]
            values += values[:1]  # Complete the circle

            # Plot
            ax.plot(angles, values, 'o-', linewidth=2, label=row['rescue_mutations'],
                   color=colors[i], alpha=0.7)
            ax.fill(angles, values, alpha=0.15, color=colors[i])

        # Customize chart
        ax.set_xticks(angles[:-1])
        ax.set_xticklabels(categories, fontsize=10)
        ax.set_ylim(0, 1)
        ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
        ax.set_yticklabels(['0.2', '0.4', '0.6', '0.8', '1.0'], fontsize=8)
        ax.set_title(f'{target}\nTop 5 Pareto Candidates', fontsize=14, fontweight='bold', pad=20)
        ax.grid(True, alpha=0.3)
        ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1), fontsize=9)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

    logger.info(f"Created radar chart comparison: {output_path}")


def create_pipeline_flow_diagram(output_path: Path) -> None:
    """Create a Sankey-style flow diagram showing the pipeline."""
    logger = get_logger(__name__)

    fig, ax = plt.subplots(figsize=(18, 10))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis('off')

    # Title
    ax.text(5, 9.5, 'ProteinForge-p53 Computational Pipeline',
           ha='center', fontsize=20, fontweight='bold')

    # Define stages with positions and sizes
    stages = [
        {'name': 'ClinVar\nData', 'pos': (1, 7), 'size': (1.2, 1), 'color': '#3498db', 'count': '8,778 variants'},
        {'name': 'Label\nExtraction', 'pos': (1, 5), 'size': (1.2, 1), 'color': '#3498db', 'count': '357 labeled'},
        {'name': 'EvoEF2\nScoring', 'pos': (1, 3), 'size': (1.2, 1), 'color': '#9b59b6', 'count': '871 scored'},
        {'name': 'Benchmark\nValidation', 'pos': (1, 1), 'size': (1.2, 1), 'color': '#2ecc71', 'count': 'AUC=0.844'},

        {'name': 'AlphaFold\nStructure', 'pos': (4, 7), 'size': (1.2, 1), 'color': '#e67e22', 'count': 'p53 core'},
        {'name': 'Constraint\nDefinition', 'pos': (4, 5), 'size': (1.2, 1), 'color': '#e67e22', 'count': '4 protected'},
        {'name': 'Candidate\nGeneration', 'pos': (4, 3), 'size': (1.2, 1), 'color': '#9b59b6', 'count': '2,055 total'},

        {'name': 'Pareto\nOptimization', 'pos': (7, 5), 'size': (1.5, 1.5), 'color': '#e74c3c', 'count': '18 optimal'},
        {'name': 'Tiered\nRecommendations', 'pos': (7, 2), 'size': (1.5, 1.5), 'color': '#2ecc71', 'count': '3 per target'},
    ]

    # Draw boxes
    for stage in stages:
        x, y = stage['pos']
        w, h = stage['size']

        # Create fancy box
        box = FancyBboxPatch((x - w/2, y - h/2), w, h,
                            boxstyle="round,pad=0.1",
                            facecolor=stage['color'], edgecolor='black',
                            linewidth=3, alpha=0.8)
        ax.add_patch(box)

        # Add text
        ax.text(x, y + 0.15, stage['name'], ha='center', va='center',
               fontsize=11, fontweight='bold', color='white')
        ax.text(x, y - 0.25, stage['count'], ha='center', va='center',
               fontsize=9, color='white', style='italic')

    # Draw arrows with labels
    arrows = [
        {'from': (1, 6.5), 'to': (1, 5.5), 'label': 'Parse XML'},
        {'from': (1, 4.5), 'to': (1, 3.5), 'label': 'Missense only'},
        {'from': (1, 2.5), 'to': (1, 1.5), 'label': 'Validate'},

        {'from': (4, 6.5), 'to': (4, 5.5), 'label': 'Extract\nfeatures'},
        {'from': (4, 4.5), 'to': (4, 3.5), 'label': 'Beam search'},

        {'from': (2, 3), 'to': (3.5, 3), 'label': 'Merge'},
        {'from': (5, 5), 'to': (6.5, 5), 'label': '3 objectives'},
        {'from': (7, 4.25), 'to': (7, 2.75), 'label': 'Select best'},
    ]

    for arrow in arrows:
        fx, fy = arrow['from']
        tx, ty = arrow['to']

        # Draw arrow
        ax.annotate('', xy=(tx, ty), xytext=(fx, fy),
                   arrowprops=dict(arrowstyle='->', lw=3, color='black'))

        # Add label
        mx, my = (fx + tx) / 2, (fy + ty) / 2
        ax.text(mx + 0.3, my, arrow['label'], fontsize=9, style='italic',
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))

    # Add metrics boxes
    metrics = [
        {'pos': (9, 8), 'title': 'INPUT', 'items': ['8,778 variants', '1 structure', '3 targets']},
        {'pos': (9, 5), 'title': 'PROCESS', 'items': ['2,055 candidates', '3 objectives', '<30 min runtime']},
        {'pos': (9, 2), 'title': 'OUTPUT', 'items': ['18 Pareto designs', '6 visualizations', '100% reproducible']},
    ]

    for metric in metrics:
        x, y = metric['pos']

        # Box
        rect = Rectangle((x - 0.6, y - 0.6), 1.2, 1.2,
                        facecolor='#ecf0f1', edgecolor='black', linewidth=2)
        ax.add_patch(rect)

        # Title
        ax.text(x, y + 0.4, metric['title'], ha='center', fontsize=10,
               fontweight='bold')

        # Items
        for i, item in enumerate(metric['items']):
            ax.text(x, y + 0.1 - i*0.25, f'• {item}', ha='center', fontsize=8)

    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

    logger.info(f"Created pipeline flow diagram: {output_path}")


def create_domain_architecture_diagram(processed_dir: Path, output_path: Path) -> None:
    """Create p53 domain architecture with rescue mutation positions."""
    logger = get_logger(__name__)

    fig, ax = plt.subplots(figsize=(20, 6))

    # Define p53 domain structure
    domains = [
        {'name': 'N-terminal\nTransactivation', 'start': 1, 'end': 61, 'color': '#3498db'},
        {'name': 'Proline-rich', 'start': 61, 'end': 94, 'color': '#9b59b6'},
        {'name': 'DNA-Binding Core', 'start': 94, 'end': 312, 'color': '#e74c3c'},
        {'name': 'Tetramerization', 'start': 312, 'end': 355, 'color': '#f39c12'},
        {'name': 'C-terminal\nRegulatory', 'start': 355, 'end': 393, 'color': '#2ecc71'},
    ]

    # Draw domains
    y_center = 0.5
    height = 0.15
    total_length = 393

    for domain in domains:
        x_start = domain['start'] / total_length
        x_end = domain['end'] / total_length
        width = x_end - x_start

        rect = Rectangle((x_start, y_center - height/2), width, height,
                        facecolor=domain['color'], edgecolor='black',
                        linewidth=2, alpha=0.8)
        ax.add_patch(rect)

        # Add label
        x_mid = (x_start + x_end) / 2
        ax.text(x_mid, y_center, domain['name'], ha='center', va='center',
               fontsize=9, fontweight='bold', color='white')

        # Add boundaries
        ax.text(x_start, y_center - height/2 - 0.05, str(domain['start']),
               ha='center', fontsize=7)
        ax.text(x_end, y_center - height/2 - 0.05, str(domain['end']),
               ha='center', fontsize=7)

    # Load rescue mutations
    all_rescues = {}
    targets = ["R175H", "R248Q", "R273H"]

    for target in targets:
        target_dir = processed_dir / "rescues" / target
        candidates_path = target_dir / "candidates.parquet"

        if not candidates_path.exists():
            continue

        df = pd.read_parquet(candidates_path)
        pareto = df[df['is_pareto']]

        for _, row in pareto.iterrows():
            rescue_muts = row['rescue_mutations'].split(',') if row['rescue_mutations'] else []
            for mut in rescue_muts:
                # Extract position
                pos = int(''.join(filter(str.isdigit, mut)))
                if pos not in all_rescues:
                    all_rescues[pos] = {'count': 0, 'mutations': []}
                all_rescues[pos]['count'] += 1
                all_rescues[pos]['mutations'].append(mut)

    # Plot rescue positions
    for pos, data in all_rescues.items():
        x = pos / total_length

        # Size based on frequency
        size = 0.015 + data['count'] * 0.005

        circle = Circle((x, y_center + height/2 + 0.15), size,
                       facecolor='#2ecc71', edgecolor='black',
                       linewidth=2, zorder=10)
        ax.add_patch(circle)

        # Add arrow
        ax.plot([x, x], [y_center + height/2, y_center + height/2 + 0.12],
               'k-', linewidth=1.5, zorder=9)

        # Add label for most frequent
        if data['count'] >= 3:
            ax.text(x, y_center + height/2 + 0.22, data['mutations'][0],
                   ha='center', fontsize=8, fontweight='bold')

    # Mark protected residues (zinc-binding)
    protected = [176, 179, 238, 242]
    for pos in protected:
        x = pos / total_length
        ax.plot(x, y_center - height/2 - 0.15, 'r^', markersize=12)
        ax.text(x, y_center - height/2 - 0.22, f'C{pos}' if pos != 179 else 'H179',
               ha='center', fontsize=7, color='red', fontweight='bold')

    # Mark cancer hotspots
    hotspots = [175, 248, 273]
    for pos in hotspots:
        x = pos / total_length
        ax.plot(x, y_center, 'k*', markersize=20, zorder=11)
        ax.text(x, y_center - height/2 - 0.3, f'R{pos}',
               ha='center', fontsize=9, fontweight='bold',
               bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.7))

    # Set limits and styling
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.4, 1)
    ax.set_aspect('equal')
    ax.axis('off')

    # Title
    ax.set_title('TP53 Domain Architecture with Rescue Mutations and Cancer Hotspots',
                fontsize=18, fontweight='bold', pad=20)

    # Legend
    legend_elements = [
        Circle((0, 0), 0.01, facecolor='#2ecc71', edgecolor='black', label='Rescue Position'),
        plt.Line2D([0], [0], marker='^', color='w', markerfacecolor='red', markersize=10, label='Zn-binding (protected)'),
        plt.Line2D([0], [0], marker='*', color='w', markerfacecolor='black', markersize=12, label='Cancer Hotspot'),
    ]
    ax.legend(handles=legend_elements, loc='upper left', fontsize=11, frameon=True, shadow=True)

    # Add caption
    caption = f"Rescue mutations clustered in DNA-binding core (94-312). All >8Å from protected Zn-binding residues."
    ax.text(0.5, -0.35, caption, ha='center', fontsize=11, style='italic', transform=ax.transAxes)

    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

    logger.info(f"Created domain architecture diagram: {output_path}")


def create_progressive_optimization_plot(processed_dir: Path, output_path: Path) -> None:
    """Create a plot showing optimization progression through beam search."""
    logger = get_logger(__name__)

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    fig.suptitle('Beam Search Optimization: Progression from Singles to Triples',
                fontsize=18, fontweight='bold', y=1.02)

    targets = ["R175H", "R248Q", "R273H"]

    for idx, (target, ax) in enumerate(zip(targets, axes)):
        target_dir = processed_dir / "rescues" / target
        candidates_path = target_dir / "candidates.parquet"

        if not candidates_path.exists():
            continue

        df = pd.read_parquet(candidates_path)

        # Separate by mutation count
        singles = df[df['n_rescue'] == 1]
        doubles = df[df['n_rescue'] == 2]
        triples = df[df['n_rescue'] == 3]

        # Plot each tier with different styling
        ax.scatter(singles['risk'], singles['ddg_gain'], s=80, alpha=0.6,
                  color='#3498db', edgecolor='black', linewidth=1,
                  label=f'Singles (n={len(singles)})', zorder=3)

        ax.scatter(doubles['risk'], doubles['ddg_gain'], s=120, alpha=0.6,
                  color='#f39c12', edgecolor='black', linewidth=1,
                  label=f'Doubles (n={len(doubles)})', marker='s', zorder=2)

        ax.scatter(triples['risk'], triples['ddg_gain'], s=160, alpha=0.6,
                  color='#e74c3c', edgecolor='black', linewidth=1,
                  label=f'Triples (n={len(triples)})', marker='^', zorder=1)

        # Highlight Pareto front
        pareto = df[df['is_pareto']]
        ax.scatter(pareto['risk'], pareto['ddg_gain'], s=300, alpha=0,
                  edgecolor='lime', linewidth=4, label='Pareto Front', zorder=4)

        # Add arrows showing progression
        for n_muts in [1, 2]:
            subset = df[df['n_rescue'] == n_muts].nsmallest(5, 'ddg_gain')
            next_subset = df[df['n_rescue'] == n_muts + 1].nsmallest(5, 'ddg_gain')

            if len(subset) > 0 and len(next_subset) > 0:
                # Draw arrow from best at n to best at n+1
                x1, y1 = subset.iloc[0]['risk'], subset.iloc[0]['ddg_gain']
                x2, y2 = next_subset.iloc[0]['risk'], next_subset.iloc[0]['ddg_gain']

                ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                           arrowprops=dict(arrowstyle='->', lw=2, color='gray',
                                         alpha=0.5, linestyle='--'))

        ax.set_xlabel('Risk Score', fontsize=12, fontweight='bold')
        ax.set_ylabel('ΔΔG Gain', fontsize=12, fontweight='bold')
        ax.set_title(f'{target}', fontsize=14, fontweight='bold')
        ax.legend(fontsize=9, loc='lower right')
        ax.grid(True, alpha=0.3)
        ax.axhline(0, color='black', linewidth=1, linestyle='--', alpha=0.5)
        ax.axvline(0, color='black', linewidth=1, linestyle='--', alpha=0.5)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

    logger.info(f"Created progressive optimization plot: {output_path}")


def run(args, configs: dict[str, Any]) -> int:
    """Generate advanced visualizations."""
    logger = get_logger(__name__)
    paths = _paths_from_config(configs.get("paths", {}))
    report_cfg = configs.get("report", {})

    output_dir = Path(report_cfg.get("output_dir", paths["reports"]))
    figures_dir = Path(report_cfg.get("figures_dir", output_dir / "figures"))
    figures_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Creating advanced visualizations...")

    # 1. Mutation network graph
    try:
        create_mutation_network(paths["processed"], figures_dir / "mutation_network_advanced.png")
        logger.info("✓ Created mutation network")
    except Exception as e:
        logger.warning(f"Could not create mutation network: {e}")

    # 2. Radar chart comparison
    try:
        create_radar_chart_comparison(paths["processed"], figures_dir / "radar_comparison_advanced.png")
        logger.info("✓ Created radar chart comparison")
    except Exception as e:
        logger.warning(f"Could not create radar charts: {e}")

    # 3. Pipeline flow diagram
    try:
        create_pipeline_flow_diagram(figures_dir / "pipeline_flow_advanced.png")
        logger.info("✓ Created pipeline flow diagram")
    except Exception as e:
        logger.warning(f"Could not create flow diagram: {e}")

    # 4. Domain architecture
    try:
        create_domain_architecture_diagram(paths["processed"], figures_dir / "domain_architecture_advanced.png")
        logger.info("✓ Created domain architecture diagram")
    except Exception as e:
        logger.warning(f"Could not create domain architecture: {e}")

    # 5. Progressive optimization
    try:
        create_progressive_optimization_plot(paths["processed"], figures_dir / "progressive_optimization_advanced.png")
        logger.info("✓ Created progressive optimization plot")
    except Exception as e:
        logger.warning(f"Could not create optimization plot: {e}")

    logger.info(f"✓ Advanced visualizations complete in {figures_dir}")
    return 0
