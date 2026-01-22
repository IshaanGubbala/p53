#!/usr/bin/env python3
"""Generate PyMOL scripts to visualize rescue mutations.

Creates PyMOL scripts that highlight the target mutation and rescue
mutations in the p53 structure.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

# Add repo root to Python path
repo_root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(repo_root))


def parse_mutation(mut_str: str) -> tuple[str, int, str]:
    """Parse mutation string like 'S95A' into (ref, pos, alt)."""
    match = re.match(r'([A-Z*])(\d+)([A-Z*])', mut_str)
    if not match:
        raise ValueError(f"Invalid mutation: {mut_str}")
    ref, pos, alt = match.groups()
    return ref, int(pos), alt


def generate_pymol_script(target: str, rescue_muts: list[str], tier: str, pdb_path: Path, output_dir: Path) -> Path:
    """Generate a PyMOL script to visualize the rescue mutations."""
    _, target_pos, _ = parse_mutation(target)

    rescue_positions = [parse_mutation(m)[1] for m in rescue_muts]

    script_name = f"visualize_{target}_{tier}.pml"
    script_path = output_dir / script_name

    # Generate PyMOL commands
    commands = [
        "# PyMOL script to visualize rescue mutations",
        f"# Target: {target}, Tier: {tier}",
        f"# Rescue mutations: {', '.join(rescue_muts)}",
        "",
        "# Load structure",
        f"load {pdb_path}",
        "hide everything",
        "",
        "# Show protein as cartoon",
        "show cartoon, all",
        "color gray80, all",
        "",
        "# Highlight target mutation in red",
        f"select target, resi {target_pos}",
        "show sticks, target",
        "color red, target",
        "label target and name CA, 'Target: {}'".format(target),
        "",
        "# Highlight rescue mutations in green",
    ]

    rescue_selection = "+".join(str(pos) for pos in rescue_positions)
    commands.extend([
        f"select rescue, resi {rescue_selection}",
        "show sticks, rescue",
        "color green, rescue",
    ])

    # Label each rescue mutation
    for i, (mut, pos) in enumerate(zip(rescue_muts, rescue_positions), 1):
        commands.append(f"label resi {pos} and name CA, 'Rescue {i}: {mut}'")

    commands.extend([
        "",
        "# Show nearby residues for context",
        f"select nearby, (resi {target_pos} around 8) or (resi {rescue_selection} around 8)",
        "show lines, nearby and not (target or rescue)",
        "",
        "# Center view on target and rescues",
        f"zoom resi {target_pos}+{rescue_selection}",
        "",
        "# Set ray tracing for high-quality image",
        "set ray_trace_mode, 1",
        "set ray_shadows, 0",
        "",
        "# Save image",
        f"png {output_dir}/structure_{target}_{tier}.png, width=1200, height=1200, dpi=300, ray=1",
        "",
        "print 'Visualization complete!'",
    ])

    script_content = "\n".join(commands)
    script_path.write_text(script_content, encoding="utf-8")

    return script_path


def main() -> None:
    reports_dir = repo_root / "reports"
    tables_dir = reports_dir / "tables"
    figures_dir = reports_dir / "figures"
    scripts_dir = reports_dir / "pymol_scripts"

    scripts_dir.mkdir(parents=True, exist_ok=True)

    # Load tiered recommendations
    recommendations_path = tables_dir / "tiered_recommendations.json"
    if not recommendations_path.exists():
        print("Error: tiered_recommendations.json not found")
        print("Run: python scripts/create_tiered_rescues.py")
        return

    recommendations = json.loads(recommendations_path.read_text(encoding="utf-8"))

    # Find PDB file
    pdb_path = repo_root / "Data" / "raw" / "alphafold"
    pdb_files = list(pdb_path.glob("*.pdb"))
    if not pdb_files:
        print(f"Error: No PDB files found in {pdb_path}")
        return

    pdb_file = pdb_files[0]  # Use first PDB file found

    print("=== Generating PyMOL Structure Visualization Scripts ===\n")
    print(f"Using structure: {pdb_file.name}\n")

    generated_scripts = []

    for target, tiers in recommendations.items():
        print(f"Processing {target}...")

        for tier_key, tier_data in tiers.items():
            rescue_muts = tier_data["rescue_mutations"]
            tier_name = tier_key.replace("best_", "")

            script_path = generate_pymol_script(
                target, rescue_muts, tier_name, pdb_file, figures_dir
            )

            generated_scripts.append(script_path)
            print(f"  {tier_key}: {', '.join(rescue_muts)}")
            print(f"    → {script_path.name}")

    print(f"\n✓ Generated {len(generated_scripts)} PyMOL scripts in {scripts_dir}")
    print("\nTo generate structure images, run PyMOL with each script:")
    print("  pymol -c script_name.pml")
    print("\nOr open PyMOL and run: @script_name.pml")

    # Create a master script to run all visualizations
    master_script = scripts_dir / "visualize_all.pml"
    master_commands = [
        "# Master script to generate all structure visualizations",
        "# Run with: pymol -c visualize_all.pml",
        "",
    ]

    for script in generated_scripts:
        master_commands.append(f"@{script}")
        master_commands.append("reinitialize")
        master_commands.append("")

    master_script.write_text("\n".join(master_commands), encoding="utf-8")
    print(f"\n✓ Master script: {master_script}")
    print(f"  Run all visualizations: pymol -c {master_script}")


if __name__ == "__main__":
    main()
