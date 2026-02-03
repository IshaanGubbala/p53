#!/usr/bin/env python3
"""Render actual 3D protein structures with rescue mutations highlighted.

Uses matplotlib for 3D rendering of PDB structures.
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import json

# Add repo root to path
repo_root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(repo_root))


def parse_pdb(pdb_path: Path) -> dict[str, Any]:
    """Parse PDB file and extract CA atom coordinates."""
    atoms = []

    with pdb_path.open('r') as f:
        for line in f:
            if line.startswith('ATOM') and ' CA ' in line:
                # Parse PDB line
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain = line[21].strip()
                res_num = int(line[22:26].strip())
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())

                atoms.append({
                    'res_num': res_num,
                    'res_name': res_name,
                    'x': x,
                    'y': y,
                    'z': z,
                    'atom': atom_name,
                    'chain': chain
                })

    return {'atoms': atoms}


def parse_mutation(mut_str: str) -> tuple[str, int, str]:
    """Parse mutation string like 'S95A' into (ref, pos, alt)."""
    import re
    match = re.match(r'([A-Z*])(\d+)([A-Z*])', mut_str)
    if not match:
        raise ValueError(f"Invalid mutation: {mut_str}")
    ref, pos, alt = match.groups()
    return ref, int(pos), alt


def render_protein_3d(pdb_data: dict, target: str, rescue_muts: list[str],
                     tier: str, output_path: Path) -> None:
    """Render 3D protein structure with mutations highlighted."""
    atoms = pdb_data['atoms']

    # Extract coordinates
    coords = np.array([[a['x'], a['y'], a['z']] for a in atoms])
    res_nums = np.array([a['res_num'] for a in atoms])

    # Parse mutations
    _, target_pos, _ = parse_mutation(target)
    rescue_positions = [parse_mutation(m)[1] for m in rescue_muts]

    # Create figure
    fig = plt.figure(figsize=(14, 12))
    ax = fig.add_subplot(111, projection='3d')

    # Plot backbone as line
    ax.plot(coords[:, 0], coords[:, 1], coords[:, 2],
           'gray', linewidth=1, alpha=0.3, zorder=1)

    # Plot all CA atoms as small dots
    ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2],
              c='lightgray', s=10, alpha=0.5, zorder=2)

    # Highlight target mutation in RED
    target_mask = res_nums == target_pos
    if target_mask.any():
        target_coords = coords[target_mask]
        ax.scatter(target_coords[:, 0], target_coords[:, 1], target_coords[:, 2],
                  c='red', s=500, alpha=0.9, edgecolor='black', linewidth=2,
                  marker='o', label=f'Target: {target}', zorder=10)

        # Add label
        ax.text(target_coords[0, 0], target_coords[0, 1], target_coords[0, 2] + 5,
               target, fontsize=12, fontweight='bold', color='red',
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    # Highlight rescue mutations in GREEN
    colors_green = plt.cm.Greens(np.linspace(0.5, 0.9, len(rescue_muts)))

    for i, (mut, color) in enumerate(zip(rescue_muts, colors_green)):
        _, pos, _ = parse_mutation(mut)
        rescue_mask = res_nums == pos

        if rescue_mask.any():
            rescue_coords = coords[rescue_mask]
            ax.scatter(rescue_coords[:, 0], rescue_coords[:, 1], rescue_coords[:, 2],
                      c=[color], s=400, alpha=0.9, edgecolor='black', linewidth=2,
                      marker='o', label=f'Rescue {i+1}: {mut}', zorder=9)

            # Add label
            label_offset = 5 + i * 3
            ax.text(rescue_coords[0, 0], rescue_coords[0, 1],
                   rescue_coords[0, 2] + label_offset,
                   mut, fontsize=10, fontweight='bold', color='green',
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    # Highlight nearby residues (within 10Å)
    all_mut_positions = [target_pos] + rescue_positions
    nearby_indices = []

    for pos in all_mut_positions:
        pos_mask = res_nums == pos
        if pos_mask.any():
            pos_coord = coords[pos_mask][0]
            distances = np.linalg.norm(coords - pos_coord, axis=1)
            nearby = distances < 10.0
            nearby_indices.extend(np.where(nearby)[0].tolist())

    nearby_indices = list(set(nearby_indices))
    if nearby_indices:
        nearby_coords = coords[nearby_indices]
        ax.scatter(nearby_coords[:, 0], nearby_coords[:, 1], nearby_coords[:, 2],
                  c='yellow', s=50, alpha=0.4, zorder=3)

    # Styling
    ax.set_xlabel('X (Å)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Y (Å)', fontsize=11, fontweight='bold')
    ax.set_zlabel('Z (Å)', fontsize=11, fontweight='bold')

    title = f'TP53 Structure: {target} + {", ".join(rescue_muts)}\n(Tier: {tier})'
    ax.set_title(title, fontsize=14, fontweight='bold', pad=20)

    # Legend
    ax.legend(loc='upper right', fontsize=9, frameon=True, shadow=True)

    # Set viewing angle
    ax.view_init(elev=20, azim=45)

    # Equal aspect ratio
    max_range = np.array([
        coords[:, 0].max() - coords[:, 0].min(),
        coords[:, 1].max() - coords[:, 1].min(),
        coords[:, 2].max() - coords[:, 2].min()
    ]).max() / 2.0

    mid_x = (coords[:, 0].max() + coords[:, 0].min()) * 0.5
    mid_y = (coords[:, 1].max() + coords[:, 1].min()) * 0.5
    mid_z = (coords[:, 2].max() + coords[:, 2].min()) * 0.5

    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)

    # Grid
    ax.grid(True, alpha=0.3)

    # Add caption
    caption = f"3D structure showing target mutation (red) and rescue mutations (green).\nYellow spheres indicate residues within 10Å."
    fig.text(0.5, 0.02, caption, ha='center', fontsize=10, style='italic', wrap=True)

    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()


def render_multiple_views(pdb_data: dict, target: str, rescue_muts: list[str],
                         tier: str, ddg_gain: float, risk: float, output_path: Path) -> None:
    """Render protein from multiple viewing angles."""
    atoms = pdb_data['atoms']
    coords = np.array([[a['x'], a['y'], a['z']] for a in atoms])
    res_nums = np.array([a['res_num'] for a in atoms])

    # Parse mutations
    _, target_pos, _ = parse_mutation(target)
    rescue_positions = [parse_mutation(m)[1] for m in rescue_muts]

    # Create figure with 4 subplots
    fig = plt.figure(figsize=(20, 16))

    views = [
        (20, 45, 'Front View'),
        (20, 135, 'Side View'),
        (80, 45, 'Top View'),
        (10, 225, 'Back View')
    ]

    for idx, (elev, azim, view_name) in enumerate(views, 1):
        ax = fig.add_subplot(2, 2, idx, projection='3d')

        # Plot backbone
        ax.plot(coords[:, 0], coords[:, 1], coords[:, 2],
               'gray', linewidth=1.5, alpha=0.4, zorder=1)

        # Plot all CA atoms
        ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2],
                  c='lightgray', s=20, alpha=0.6, zorder=2)

        # Highlight target
        target_mask = res_nums == target_pos
        if target_mask.any():
            target_coords = coords[target_mask]
            ax.scatter(target_coords[:, 0], target_coords[:, 1], target_coords[:, 2],
                      c='red', s=600, alpha=0.95, edgecolor='black', linewidth=3,
                      marker='o', zorder=10)

        # Highlight rescues
        for mut in rescue_muts:
            _, pos, _ = parse_mutation(mut)
            rescue_mask = res_nums == pos
            if rescue_mask.any():
                rescue_coords = coords[rescue_mask]
                ax.scatter(rescue_coords[:, 0], rescue_coords[:, 1], rescue_coords[:, 2],
                          c='green', s=500, alpha=0.95, edgecolor='black', linewidth=2,
                          marker='o', zorder=9)

        # Styling
        ax.set_xlabel('X (Å)', fontsize=10)
        ax.set_ylabel('Y (Å)', fontsize=10)
        ax.set_zlabel('Z (Å)', fontsize=10)
        ax.set_title(view_name, fontsize=13, fontweight='bold')
        ax.view_init(elev=elev, azim=azim)
        ax.grid(True, alpha=0.3)

        # Equal aspect
        max_range = np.array([
            coords[:, 0].max() - coords[:, 0].min(),
            coords[:, 1].max() - coords[:, 1].min(),
            coords[:, 2].max() - coords[:, 2].min()
        ]).max() / 2.0

        mid_x = (coords[:, 0].max() + coords[:, 0].min()) * 0.5
        mid_y = (coords[:, 1].max() + coords[:, 1].min()) * 0.5
        mid_z = (coords[:, 2].max() + coords[:, 2].min()) * 0.5

        ax.set_xlim(mid_x - max_range, mid_x + max_range)
        ax.set_ylim(mid_y - max_range, mid_y + max_range)
        ax.set_zlim(mid_z - max_range, mid_z + max_range)

    # Overall title
    title = f'TP53 Structure - Multiple Views: {target} + {", ".join(rescue_muts)}'
    fig.suptitle(title, fontsize=16, fontweight='bold', y=0.98)

    # Metrics box
    metrics_text = f'Tier: {tier} | ΔΔG Gain: {ddg_gain:.1f} | Risk: {risk:.3f}'
    fig.text(0.5, 0.01, metrics_text, ha='center', fontsize=11, fontweight='bold',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    plt.tight_layout(rect=[0, 0.02, 1, 0.97])
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()


def main():
    """Generate 3D protein structure renders for all rescue designs."""
    data_dir = repo_root / "Data"
    reports_dir = repo_root / "reports"
    figures_dir = reports_dir / "figures"
    structures_dir = figures_dir / "structures_3d"

    structures_dir.mkdir(parents=True, exist_ok=True)

    # Find PDB file
    pdb_path = data_dir / "raw" / "alphafold"
    pdb_files = list(pdb_path.glob("*.pdb"))
    if not pdb_files:
        print(f"Error: No PDB files found in {pdb_path}")
        return

    pdb_file = pdb_files[0]

    print("=" * 70)
    print("3D Protein Structure Renderer")
    print("=" * 70)
    print(f"\nParsing PDB: {pdb_file.name}")

    # Parse PDB
    pdb_data = parse_pdb(pdb_file)
    print(f"Loaded {len(pdb_data['atoms'])} CA atoms\n")

    # Load recommendations
    recommendations_path = reports_dir / "tables" / "tiered_recommendations.json"
    if not recommendations_path.exists():
        print(f"Error: {recommendations_path} not found")
        return

    recommendations = json.loads(recommendations_path.read_text())

    print("Rendering 3D structures:\n")

    # Render each design
    for target, tiers in recommendations.items():
        print(f"  {target}:")

        for tier_key, tier_data in tiers.items():
            rescue_muts = tier_data["rescue_mutations"]
            tier_name = tier_key.replace("best_", "")
            ddg_gain = tier_data["ddg_gain"]
            risk = tier_data["risk"]

            # Single view
            output_single = structures_dir / f"protein_3d_{target}_{tier_name}.png"
            render_protein_3d(pdb_data, target, rescue_muts, tier_name, output_single)
            print(f"    ✓ {output_single.name}")

            # Multiple views
            output_multi = structures_dir / f"protein_3d_{target}_{tier_name}_multiview.png"
            render_multiple_views(pdb_data, target, rescue_muts, tier_name,
                                ddg_gain, risk, output_multi)
            print(f"    ✓ {output_multi.name}")

    print("\n" + "=" * 70)
    print(f"✓ 3D structures saved to: {structures_dir}")
    print("=" * 70)
    print("\nView the rendered structures:")
    print(f"  open {structures_dir}/protein_3d_*.png")


if __name__ == "__main__":
    main()
