#!/usr/bin/env python3
"""Generate high-quality protein structure renders using PyMOL.

Uses PyMOL for ray-traced, publication-quality renders of p53 structures
with rescue mutations highlighted.
"""
from __future__ import annotations

import sys
import json
from pathlib import Path
from typing import Any

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


def generate_pymol_script(
    pdb_path: Path,
    target: str,
    rescue_muts: list[str],
    tier: str,
    output_path: Path,
    script_path: Path,
    ddg_gain: float,
    risk: float
) -> None:
    """Generate a PyMOL script for rendering a specific rescue design."""

    # Parse mutation positions
    _, target_pos, _ = parse_mutation(target)
    rescue_positions = [parse_mutation(m)[1] for m in rescue_muts]

    # Create PyMOL commands
    commands = [
        "# PyMOL script for rendering p53 rescue mutations",
        f"# Target: {target}, Tier: {tier}",
        f"# ΔΔG Gain: {ddg_gain:.1f}, Risk: {risk:.3f}",
        "",
        "# Load structure",
        f"load {pdb_path}",
        "remove solvent",
        "",
        "# Set background and basic display",
        "bg_color white",
        "hide everything",
        "show cartoon",
        "color gray80, all",
        "set cartoon_smooth_loops, 0",
        "set cartoon_flat_sheets, 0",
        "set cartoon_fancy_helices, 1",
        "set ray_shadows, 1",
        "set ray_trace_mode, 1",
        "set antialias, 2",
        "set ambient, 0.4",
        "set specular, 0.5",
        "set shininess, 10",
        "set depth_cue, 1",
        "set ray_trace_fog, 1",
        "",
        "# Highlight target mutation in RED",
        f"select target_mut, resi {target_pos}",
        "show spheres, target_mut",
        "color red, target_mut",
        "set sphere_scale, 1.5, target_mut",
        "set sphere_transparency, 0.0, target_mut",
        "",
    ]

    # Add rescue mutations
    for i, (mut, pos) in enumerate(zip(rescue_muts, rescue_positions), 1):
        commands.extend([
            f"# Rescue mutation {i}: {mut}",
            f"select rescue_{i}, resi {pos}",
            f"show spheres, rescue_{i}",
            f"color green, rescue_{i}",
            f"set sphere_scale, 1.3, rescue_{i}",
            f"set sphere_transparency, 0.0, rescue_{i}",
            "",
        ])

    # Add nearby residues
    all_positions = [target_pos] + rescue_positions
    positions_str = "+".join(str(p) for p in all_positions)

    commands.extend([
        "# Show nearby residues (within 10Å)",
        f"select mutation_sites, resi {positions_str}",
        "select nearby, (mutation_sites around 10) and not mutation_sites",
        "show sticks, nearby",
        "color yellow, nearby",
        "set stick_radius, 0.15, nearby",
        "set stick_transparency, 0.5, nearby",
        "",
        "# Set view",
        "orient",
        "zoom all, 5",
        "",
        "# Add labels",
        f"label target_mut and name CA, '{target}'",
    ])

    for i, mut in enumerate(rescue_muts, 1):
        commands.append(f"label rescue_{i} and name CA, '{mut}'")

    commands.extend([
        "",
        "set label_size, 20",
        "set label_color, black",
        "set label_bg_color, white",
        "set label_bg_transparency, 0.3",
        "",
        "# Ray trace and save",
        "viewport 2400, 2000",
        f"ray 2400, 2000",
        f"png {output_path}, dpi=300",
        "",
        "# Deselect all",
        "deselect",
        "",
        "print('Render complete!')",
    ])

    # Write script
    script_path.write_text("\n".join(commands))


def generate_multiview_script(
    pdb_path: Path,
    target: str,
    rescue_muts: list[str],
    tier: str,
    output_dir: Path,
    script_path: Path,
    ddg_gain: float,
    risk: float
) -> list[Path]:
    """Generate PyMOL script for multiple viewing angles."""

    # Parse mutation positions
    _, target_pos, _ = parse_mutation(target)
    rescue_positions = [parse_mutation(m)[1] for m in rescue_muts]

    # Define views
    views = [
        ('front', 'Front View'),
        ('side', 'Side View'),
        ('top', 'Top View'),
        ('back', 'Back View'),
    ]

    output_paths = []

    # Create PyMOL commands
    commands = [
        "# PyMOL script for multi-view rendering",
        f"# Target: {target}, Tier: {tier}",
        "",
        "# Load structure",
        f"load {pdb_path}",
        "remove solvent",
        "",
        "# Set background and basic display",
        "bg_color white",
        "hide everything",
        "show cartoon",
        "color gray80, all",
        "set cartoon_smooth_loops, 0",
        "set cartoon_fancy_helices, 1",
        "set ray_shadows, 1",
        "set antialias, 2",
        "set ambient, 0.4",
        "set specular, 0.5",
        "",
        "# Highlight target mutation",
        f"select target_mut, resi {target_pos}",
        "show spheres, target_mut",
        "color red, target_mut",
        "set sphere_scale, 1.5, target_mut",
        "",
    ]

    # Add rescue mutations
    for i, (mut, pos) in enumerate(zip(rescue_muts, rescue_positions), 1):
        commands.extend([
            f"select rescue_{i}, resi {pos}",
            f"show spheres, rescue_{i}",
            f"color green, rescue_{i}",
            f"set sphere_scale, 1.3, rescue_{i}",
        ])

    # Add nearby residues
    all_positions = [target_pos] + rescue_positions
    positions_str = "+".join(str(p) for p in all_positions)

    commands.extend([
        f"select mutation_sites, resi {positions_str}",
        "select nearby, (mutation_sites around 10) and not mutation_sites",
        "show sticks, nearby",
        "color yellow, nearby",
        "set stick_radius, 0.15, nearby",
        "set stick_transparency, 0.5, nearby",
        "",
        "# Orient",
        "orient",
        "zoom all, 5",
        "",
        "viewport 1200, 1000",
        "",
    ])

    # Add rendering commands for each view
    for view_name, view_label in views:
        output_path = output_dir / f"protein_pymol_{target}_{tier}_{view_name}.png"
        output_paths.append(output_path)

        commands.extend([
            f"# {view_label}",
        ])

        # Different rotations
        if view_name == 'side':
            commands.append("turn y, 90")
        elif view_name == 'top':
            commands.append("turn x, 90")
        elif view_name == 'back':
            commands.append("turn y, 180")

        commands.extend([
            f"ray 1200, 1000",
            f"png {output_path}, dpi=300",
            "",
        ])

        # Reset view
        if view_name != 'back':
            commands.extend([
                "orient",
                "zoom all, 5",
                "",
            ])

    commands.append("print('Multi-view render complete!')")

    # Write script
    script_path.write_text("\n".join(commands))
    return output_paths


def main():
    """Generate PyMOL renders for all rescue designs."""
    data_dir = repo_root / "Data"
    reports_dir = repo_root / "reports"
    figures_dir = reports_dir / "figures"
    pymol_dir = figures_dir / "pymol_renders"
    scripts_dir = reports_dir / "pymol_scripts"

    pymol_dir.mkdir(parents=True, exist_ok=True)
    scripts_dir.mkdir(parents=True, exist_ok=True)

    # Find PDB file
    pdb_path = data_dir / "raw" / "alphafold"
    pdb_files = list(pdb_path.glob("*.pdb"))
    if not pdb_files:
        print(f"Error: No PDB files found in {pdb_path}")
        return

    pdb_file = pdb_files[0]

    print("=" * 70)
    print("PyMOL High-Quality Protein Structure Renderer")
    print("=" * 70)
    print(f"\nUsing PDB: {pdb_file.name}")
    print(f"PyMOL path: /Applications/PyMOL.app/Contents/MacOS/PyMOL\n")

    # Load recommendations
    recommendations_path = reports_dir / "tables" / "tiered_recommendations.json"
    if not recommendations_path.exists():
        print(f"Error: {recommendations_path} not found")
        return

    recommendations = json.loads(recommendations_path.read_text())

    print("Generating PyMOL scripts:\n")

    pymol_exe = "/Applications/PyMOL.app/Contents/MacOS/PyMOL"
    all_scripts = []

    # Generate scripts for each design
    for target, tiers in recommendations.items():
        print(f"  {target}:")

        for tier_key, tier_data in tiers.items():
            rescue_muts = tier_data["rescue_mutations"]
            tier_name = tier_key.replace("best_", "")
            ddg_gain = tier_data["ddg_gain"]
            risk = tier_data["risk"]

            # Single high-quality render
            output_single = pymol_dir / f"protein_pymol_{target}_{tier_name}_hq.png"
            script_single = scripts_dir / f"render_{target}_{tier_name}_hq.pml"

            generate_pymol_script(
                pdb_file, target, rescue_muts, tier_name,
                output_single, script_single, ddg_gain, risk
            )
            all_scripts.append(script_single)
            print(f"    ✓ Script: {script_single.name}")

            # Multi-view renders
            script_multi = scripts_dir / f"render_{target}_{tier_name}_multiview.pml"
            output_paths = generate_multiview_script(
                pdb_file, target, rescue_muts, tier_name,
                pymol_dir, script_multi, ddg_gain, risk
            )
            all_scripts.append(script_multi)
            print(f"    ✓ Script: {script_multi.name}")

    print(f"\n{'=' * 70}")
    print(f"Generated {len(all_scripts)} PyMOL scripts")
    print(f"{'=' * 70}\n")

    print("Executing PyMOL scripts (this may take several minutes)...\n")

    # Execute each script
    for i, script in enumerate(all_scripts, 1):
        print(f"  [{i}/{len(all_scripts)}] Rendering {script.name}...")

        import subprocess
        result = subprocess.run(
            [pymol_exe, "-cq", str(script)],
            capture_output=True,
            text=True
        )

        if result.returncode == 0:
            print(f"    ✓ Success")
        else:
            print(f"    ✗ Error: {result.stderr[:100]}")

    print(f"\n{'=' * 70}")
    print(f"✓ All renders saved to: {pymol_dir}")
    print(f"{'=' * 70}\n")

    # Count output files
    png_files = list(pymol_dir.glob("protein_pymol_*.png"))
    print(f"Generated {len(png_files)} high-quality images:\n")

    for target in ["R175H", "R248Q", "R273H"]:
        target_files = [f for f in png_files if target in f.name]
        if target_files:
            print(f"  {target}:")
            for f in sorted(target_files):
                print(f"    - {f.name}")

    print(f"\nView renders:")
    print(f"  open {pymol_dir}/protein_pymol_*.png")


if __name__ == "__main__":
    main()
