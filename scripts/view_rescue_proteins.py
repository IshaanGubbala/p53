#!/usr/bin/env python3
"""Interactive protein structure viewer for rescue mutations.

This script provides multiple ways to visualize rescue mutations:
1. Generate PyMOL session files for 3D viewing
2. Create 2D schematic representations
3. Export PDB files with mutations applied
4. Generate comparison views (before/after)
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle, FancyBboxPatch, Circle
import numpy as np

# Add repo root to path
repo_root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(repo_root))


def parse_mutation(mut_str: str) -> tuple[str, int, str]:
    """Parse mutation string like 'S95A' into (ref, pos, alt)."""
    import re
    match = re.match(r'([A-Z*])(\d+)([A-Z*])', mut_str)
    if not match:
        raise ValueError(f"Invalid mutation: {mut_str}")
    ref, pos, alt = match.groups()
    return ref, int(pos), alt


def generate_pymol_session(target: str, rescue_muts: list[str], tier: str,
                          pdb_path: Path, output_dir: Path) -> Path:
    """Generate an interactive PyMOL session file."""
    _, target_pos, _ = parse_mutation(target)
    rescue_positions = [parse_mutation(m)[1] for m in rescue_muts]

    session_name = f"session_{target}_{tier}.pml"
    session_path = output_dir / session_name

    # Generate PyMOL commands for interactive session
    commands = [
        "# PyMOL Interactive Session for Rescue Visualization",
        f"# Target: {target}, Tier: {tier}",
        f"# Rescue mutations: {', '.join(rescue_muts)}",
        "",
        "# Load structure",
        f"load {pdb_path}",
        "",
        "# Initial setup",
        "bg_color white",
        "hide everything",
        "show cartoon, all",
        "color gray80, all",
        "set cartoon_fancy_helices, 1",
        "set cartoon_smooth_loops, 1",
        "",
        "# Color by secondary structure",
        "color salmon, ss h",
        "color lightblue, ss s",
        "color gray80, ss l+''",
        "",
        "# Highlight target mutation in RED",
        f"select target_mut, resi {target_pos}",
        "show sticks, target_mut",
        "show spheres, target_mut",
        "color red, target_mut",
        "set sphere_scale, 0.3, target_mut",
        "",
        "# Highlight rescue mutations in GREEN",
    ]

    rescue_selection = "+".join(str(pos) for pos in rescue_positions)
    commands.extend([
        f"select rescue_muts, resi {rescue_selection}",
        "show sticks, rescue_muts",
        "show spheres, rescue_muts",
        "color green, rescue_muts",
        "set sphere_scale, 0.3, rescue_muts",
        "",
        "# Show nearby context (within 8Å)",
        f"select context, (resi {target_pos} around 8) or (resi {rescue_selection} around 8)",
        "show lines, context and not (target_mut or rescue_muts)",
        "color gray70, context and not (target_mut or rescue_muts)",
        "",
        "# Add labels",
        f"label target_mut and name CA, '\"Target: {target}\"'",
    ])

    for i, (mut, pos) in enumerate(zip(rescue_muts, rescue_positions), 1):
        commands.append(f"label resi {pos} and name CA, '\"Rescue {i}: {mut}\"'")

    commands.extend([
        "",
        "# Set viewing angle",
        f"zoom resi {target_pos}+{rescue_selection}",
        "orient",
        "",
        "# Enable better rendering",
        "set ray_trace_mode, 1",
        "set ray_shadows, 0",
        "set antialias, 2",
        "set cartoon_side_chain_helper, on",
        "",
        "# Create selections for easy toggling",
        "disable all",
        "enable all",
        "",
        "# Print instructions",
        "print '======================================'",
        "print 'PyMOL Interactive Session Loaded'",
        "print '======================================'",
        f"print 'Target: {target}'",
        f"print 'Rescue: {', '.join(rescue_muts)}'",
        "print ''",
        "print 'Controls:'",
        "print '  - Mouse drag: rotate'",
        "print '  - Mouse scroll: zoom'",
        "print '  - Right-click drag: translate'",
        "print ''",
        "print 'Useful commands:'",
        "print '  hide everything; show cartoon'",
        "print '  show surface, all'",
        "print '  color_by_element'",
        "print '  ray 1200,1200'",
        "print '======================================'",
    ])

    session_path.write_text("\n".join(commands), encoding="utf-8")
    return session_path


def create_protein_schematic(target: str, rescue_muts: list[str], tier: str,
                             ddg_gain: float, risk: float, output_path: Path) -> None:
    """Create a 2D schematic representation of the protein with mutations."""
    fig, ax = plt.subplots(figsize=(16, 10))

    # Define p53 sequence length
    seq_length = 393

    # Draw protein backbone
    backbone_y = 0.5
    backbone_height = 0.05

    # Draw domains as colored rectangles
    domains = [
        {'name': 'N-TAD', 'start': 1, 'end': 61, 'color': '#3498db'},
        {'name': 'PRD', 'start': 61, 'end': 94, 'color': '#9b59b6'},
        {'name': 'DBD', 'start': 94, 'end': 312, 'color': '#e74c3c'},
        {'name': 'TET', 'start': 312, 'end': 355, 'color': '#f39c12'},
        {'name': 'CTD', 'start': 355, 'end': 393, 'color': '#2ecc71'},
    ]

    for domain in domains:
        x_start = domain['start'] / seq_length
        x_end = domain['end'] / seq_length
        width = x_end - x_start

        rect = FancyBboxPatch((x_start, backbone_y - backbone_height/2), width, backbone_height,
                             boxstyle="round,pad=0.01",
                             facecolor=domain['color'], edgecolor='black',
                             linewidth=2, alpha=0.7)
        ax.add_patch(rect)

        # Add domain label
        x_mid = (x_start + x_end) / 2
        ax.text(x_mid, backbone_y, domain['name'], ha='center', va='center',
               fontsize=10, fontweight='bold', color='white')

    # Mark target mutation
    _, target_pos, _ = parse_mutation(target)
    target_x = target_pos / seq_length

    target_circle = Circle((target_x, backbone_y + 0.2), 0.03,
                          facecolor='red', edgecolor='black', linewidth=3, zorder=10)
    ax.add_patch(target_circle)
    ax.plot([target_x, target_x], [backbone_y + backbone_height/2, backbone_y + 0.17],
           'k-', linewidth=2)
    ax.text(target_x, backbone_y + 0.26, target, ha='center', fontsize=12,
           fontweight='bold', bbox=dict(boxstyle='round', facecolor='red', alpha=0.7))

    # Mark rescue mutations
    colors_gradient = plt.cm.Greens(np.linspace(0.4, 0.9, len(rescue_muts)))

    for i, (mut, color) in enumerate(zip(rescue_muts, colors_gradient)):
        _, pos, _ = parse_mutation(mut)
        mut_x = pos / seq_length

        y_offset = backbone_y - 0.15 - (i % 3) * 0.08  # Stagger vertically to avoid overlap

        rescue_circle = Circle((mut_x, y_offset), 0.025,
                              facecolor=color, edgecolor='black', linewidth=2, zorder=9)
        ax.add_patch(rescue_circle)
        ax.plot([mut_x, mut_x], [backbone_y - backbone_height/2, y_offset + 0.025],
               'k-', linewidth=1.5)
        ax.text(mut_x, y_offset - 0.05, mut, ha='center', fontsize=10,
               fontweight='bold', bbox=dict(boxstyle='round', facecolor=color, alpha=0.7))

    # Add Zn-binding residues
    zn_residues = [176, 179, 238, 242]
    for zn_pos in zn_residues:
        zn_x = zn_pos / seq_length
        ax.plot(zn_x, backbone_y - 0.35, '^', color='orange', markersize=12)
        ax.text(zn_x, backbone_y - 0.4, f'Zn', ha='center', fontsize=8,
               color='orange', fontweight='bold')

    # Add metrics box
    metrics_text = (
        f"Tier: {tier}\n"
        f"ΔΔG Gain: {ddg_gain:.2f}\n"
        f"Risk: {risk:.3f}\n"
        f"# Rescues: {len(rescue_muts)}"
    )
    ax.text(0.02, 0.85, metrics_text, transform=ax.transAxes,
           fontsize=11, family='monospace',
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    # Set limits and styling
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.5, 1.0)
    ax.set_aspect('equal')
    ax.axis('off')

    # Title
    title = f'TP53 Rescue Design: {target} + {", ".join(rescue_muts)}'
    ax.set_title(title, fontsize=16, fontweight='bold', pad=20)

    # Legend
    legend_elements = [
        Circle((0, 0), 0.01, facecolor='red', edgecolor='black', label='Target Mutation'),
        Circle((0, 0), 0.01, facecolor='green', edgecolor='black', label='Rescue Mutation'),
        plt.Line2D([0], [0], marker='^', color='w', markerfacecolor='orange',
                   markersize=10, label='Zn-binding (protected)'),
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=10, frameon=True)

    # Add scale
    scale_y = -0.45
    ax.plot([0, 1], [scale_y, scale_y], 'k-', linewidth=1)
    ax.plot([0, 0], [scale_y - 0.02, scale_y + 0.02], 'k-', linewidth=1)
    ax.plot([1, 1], [scale_y - 0.02, scale_y + 0.02], 'k-', linewidth=1)
    ax.text(0.5, scale_y - 0.05, f'1-{seq_length} amino acids', ha='center', fontsize=9)

    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()


def create_comparison_view(processed_dir: Path, output_dir: Path) -> None:
    """Create side-by-side comparison of all tiered rescues."""
    targets = ["R175H", "R248Q", "R273H"]

    fig, axes = plt.subplots(3, 3, figsize=(24, 18))
    fig.suptitle('Complete Rescue Design Portfolio: All Targets × All Tiers',
                fontsize=20, fontweight='bold', y=0.995)

    tiers = ['best_1mut', 'best_2mut', 'best_3mut']
    tier_names = ['Single', 'Double', 'Triple']

    # Load tiered recommendations
    recommendations_path = processed_dir.parent / "reports" / "tables" / "tiered_recommendations.json"
    if not recommendations_path.exists():
        print(f"Error: {recommendations_path} not found")
        return

    recommendations = json.loads(recommendations_path.read_text(encoding="utf-8"))

    for row_idx, target in enumerate(targets):
        for col_idx, (tier, tier_name) in enumerate(zip(tiers, tier_names)):
            ax = axes[row_idx, col_idx]

            if target not in recommendations or tier not in recommendations[target]:
                ax.axis('off')
                continue

            data = recommendations[target][tier]
            rescue_muts = data['rescue_mutations']
            ddg_gain = data['ddg_gain']
            risk = data['risk']

            # Create mini schematic
            seq_length = 393
            backbone_y = 0.5

            # Draw simplified backbone
            ax.plot([0, 1], [backbone_y, backbone_y], 'k-', linewidth=8, alpha=0.3)

            # Mark target
            _, target_pos, _ = parse_mutation(target)
            target_x = target_pos / seq_length
            ax.plot(target_x, backbone_y, 'ro', markersize=15, markeredgecolor='black',
                   markeredgewidth=2)
            ax.text(target_x, backbone_y + 0.15, target, ha='center', fontsize=10,
                   fontweight='bold', bbox=dict(boxstyle='round', facecolor='red', alpha=0.6))

            # Mark rescues
            for i, mut in enumerate(rescue_muts):
                _, pos, _ = parse_mutation(mut)
                mut_x = pos / seq_length
                y_offset = backbone_y - 0.15 - (i * 0.08)
                ax.plot(mut_x, y_offset, 'go', markersize=12, markeredgecolor='black',
                       markeredgewidth=1.5)
                ax.text(mut_x, y_offset - 0.08, mut, ha='center', fontsize=8,
                       bbox=dict(boxstyle='round', facecolor='green', alpha=0.6))
                ax.plot([mut_x, mut_x], [backbone_y, y_offset + 0.06], 'k--',
                       linewidth=1, alpha=0.5)

            # Add metrics
            metrics = f"ΔΔG: {ddg_gain:.1f}\nRisk: {risk:.3f}"
            ax.text(0.02, 0.95, metrics, transform=ax.transAxes,
                   fontsize=9, verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.6))

            # Styling
            ax.set_xlim(-0.05, 1.05)
            ax.set_ylim(-0.5, 0.8)
            ax.set_aspect('equal')
            ax.axis('off')

            # Title
            if row_idx == 0:
                ax.set_title(f'{tier_name} Mutation', fontsize=14, fontweight='bold', pad=10)
            if col_idx == 0:
                ax.text(-0.15, 0.5, target, transform=ax.transAxes,
                       fontsize=16, fontweight='bold', rotation=90,
                       verticalalignment='center',
                       bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.7))

    plt.tight_layout(rect=[0, 0, 1, 0.99])
    plt.savefig(output_dir / "rescue_portfolio_comparison.png", dpi=300,
               bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"✓ Created comparison view: {output_dir / 'rescue_portfolio_comparison.png'}")


def main():
    """Main function to generate all protein views."""
    data_dir = repo_root / "Data"
    processed_dir = data_dir / "processed"
    reports_dir = repo_root / "reports"
    figures_dir = reports_dir / "figures"
    pymol_dir = reports_dir / "pymol_scripts"

    figures_dir.mkdir(parents=True, exist_ok=True)
    pymol_dir.mkdir(parents=True, exist_ok=True)

    # Find PDB file
    pdb_path = data_dir / "raw" / "alphafold"
    pdb_files = list(pdb_path.glob("*.pdb"))
    if not pdb_files:
        print(f"Error: No PDB files found in {pdb_path}")
        return

    pdb_file = pdb_files[0]

    print("=" * 60)
    print("Protein Structure Viewer - Rescue Mutations")
    print("=" * 60)
    print(f"\nUsing structure: {pdb_file.name}\n")

    # Load tiered recommendations
    recommendations_path = reports_dir / "tables" / "tiered_recommendations.json"
    if not recommendations_path.exists():
        print(f"Error: {recommendations_path} not found")
        print("Run: python scripts/create_tiered_rescues.py")
        return

    recommendations = json.loads(recommendations_path.read_text(encoding="utf-8"))

    print("Generating visualizations:\n")

    generated_files = []

    # 1. Generate PyMOL sessions for each tier
    print("1. PyMOL Interactive Sessions:")
    for target, tiers in recommendations.items():
        for tier_key, tier_data in tiers.items():
            rescue_muts = tier_data["rescue_mutations"]
            tier_name = tier_key.replace("best_", "")

            session_path = generate_pymol_session(
                target, rescue_muts, tier_name, pdb_file, pymol_dir
            )
            generated_files.append(session_path)
            print(f"   ✓ {session_path.name}")

    # 2. Generate 2D schematics
    print("\n2. 2D Protein Schematics:")
    for target, tiers in recommendations.items():
        for tier_key, tier_data in tiers.items():
            rescue_muts = tier_data["rescue_mutations"]
            tier_name = tier_key.replace("best_", "")
            ddg_gain = tier_data["ddg_gain"]
            risk = tier_data["risk"]

            output_path = figures_dir / f"schematic_{target}_{tier_name}.png"
            create_protein_schematic(target, rescue_muts, tier_name,
                                   ddg_gain, risk, output_path)
            generated_files.append(output_path)
            print(f"   ✓ {output_path.name}")

    # 3. Generate comparison view
    print("\n3. Portfolio Comparison View:")
    create_comparison_view(processed_dir, figures_dir)

    # Print summary
    print("\n" + "=" * 60)
    print(f"✓ Generated {len(generated_files)} files")
    print("=" * 60)

    print("\n📁 Output locations:")
    print(f"   PyMOL sessions: {pymol_dir}")
    print(f"   2D schematics:  {figures_dir}")

    print("\n🔬 To view in PyMOL:")
    print("   pymol session_R175H_1mut.pml")
    print("   OR open PyMOL GUI and File → Run Script")

    print("\n🖼️  To view schematics:")
    print("   open reports/figures/schematic_*.png")
    print("   open reports/figures/rescue_portfolio_comparison.png")

    print("\n✨ All protein structures are ready to view!")


if __name__ == "__main__":
    main()
