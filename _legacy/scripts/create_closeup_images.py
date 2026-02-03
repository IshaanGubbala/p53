#!/usr/bin/env python3
"""Create individual close-up images showing DNA binding interactions.

For each target, creates 3 separate images:
1. Wild-type - normal DNA binding
2. Cancer mutant - disrupted DNA binding
3. Rescued - restored DNA binding

All zoomed in on the mutation site and DNA interface.
"""
from __future__ import annotations

import sys
from pathlib import Path

repo_root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(repo_root))


def generate_closeup_script(
    protein_pdb: Path,
    dna_pdb: Path,
    target: str,
    state: str,  # "wildtype", "cancer", or "rescued"
    output_path: Path,
    script_path: Path,
) -> None:
    """Generate PyMOL script for close-up DNA binding image."""

    import re
    match = re.match(r'([A-Z*])(\d+)([A-Z*])', target)
    if not match:
        raise ValueError(f"Invalid target: {target}")

    _, target_pos, _ = match.groups()
    target_pos = int(target_pos)

    # DNA binding interface residues
    interface_residues = f"{target_pos}+241+248+273+280"

    # Color scheme based on state
    if state == "wildtype":
        mutation_color = "cyan"
        surface_color = "palegreen"
        hbond_color = "lime"
        hbond_width = 4
        hbond_gap = 0.1
        label = "WILD-TYPE\\nNormal DNA Binding"
        label_color = "blue"
    elif state == "cancer":
        mutation_color = "red"
        surface_color = "salmon"
        hbond_color = "red"
        hbond_width = 2
        hbond_gap = 0.4
        label = f"CANCER MUTANT\\n{target} - Disrupted Binding"
        label_color = "red"
    else:  # rescued
        mutation_color = "yellow"
        surface_color = "lightgreen"
        hbond_color = "lime"
        hbond_width = 4
        hbond_gap = 0.1
        label = f"RESCUED\\nBinding Restored"
        label_color = "green"

    commands = [
        f"# Close-up: {target} {state.title()} DNA Binding",
        "",
        "# High-quality settings",
        "set ray_trace_mode, 1",
        "set ray_shadows, 1",
        "set antialias, 2",
        "set ambient, 0.5",
        "set specular, 0.6",
        "set depth_cue, 1",
        "set bg_rgb, [1, 1, 1]",
        "set cartoon_fancy_helices, 1",
        "",
        "# Load structures",
        f"load {protein_pdb}, protein",
        f"load {dna_pdb}, dna",
        "",
        "# Align",
        "align dna, protein",
        "",
        "# Hide all initially",
        "hide everything",
        "",
        "# Show protein cartoon",
        "show cartoon, protein",
        "color gray70, protein",
        "",
        "# Show DNA as cartoon",
        "show cartoon, dna",
        "color marine, dna",
        "set cartoon_ring_mode, 3, dna",
        "set cartoon_ladder_mode, 1, dna",
        "",
        "# Highlight mutation site",
        f"select mut_site, protein and resi {target_pos}",
        "show spheres, mut_site",
        f"color {mutation_color}, mut_site",
        "set sphere_scale, 1.8, mut_site",
        "",
        "# Show DNA binding interface as sticks",
        f"select interface, protein and resi {interface_residues}",
        "show sticks, interface",
        "color orange, interface",
        "set stick_radius, 0.35, interface",
        "",
        "# Show hydrogen bonds to DNA",
        "distance hbonds, interface, dna, 3.5",
        "hide labels, hbonds",
        f"color {hbond_color}, hbonds",
        f"set dash_width, {hbond_width}, hbonds",
        f"set dash_gap, {hbond_gap}, hbonds",
        "",
        "# Show binding pocket surface",
        f"select pocket, protein and (interface around 10)",
        "show surface, pocket",
        "set transparency, 0.6, pocket",
        f"color {surface_color}, pocket",
        "",
        "# Highlight DNA bases near binding site",
        "select dna_contact, dna and (interface around 5)",
        "show sticks, dna_contact",
        "color yellow, dna_contact and elem C",
        "set stick_radius, 0.25, dna_contact",
        "",
        "# Label",
        "set label_size, 60",
        f"set label_color, {label_color}",
        "set label_bg_color, white",
        "set label_bg_transparency, 0.3",
        "pseudoatom label_obj, pos=[0,0,0]",
        f"label label_obj, '{label}'",
        "",
        "# ZOOM IN on binding interface (key difference!)",
        "center interface",
        "orient interface",
        "zoom interface, 4",
        "",
        "# Rotate to show binding clearly",
        "turn y, 20",
        "turn x, 10",
        "",
        "# Render high-quality image",
        "viewport 1920, 1080",
        "ray 1920, 1080",
        f"png {output_path}, dpi=300",
        "",
        f"print('Image saved: {output_path}')",
    ]

    script_path.write_text("\n".join(commands))


def main():
    """Generate close-up DNA binding images."""

    import json

    data_dir = repo_root / "Data"
    reports_dir = repo_root / "reports"
    figures_dir = reports_dir / "figures"
    closeups_dir = figures_dir / "dna_binding_closeups"
    scripts_dir = reports_dir / "pymol_scripts"

    closeups_dir.mkdir(parents=True, exist_ok=True)
    scripts_dir.mkdir(parents=True, exist_ok=True)

    # DNA structure
    dna_pdb = data_dir / "raw" / "experimental_pdbs" / "2OCJ_full.pdb"

    # Wild-type structure
    wt_pdb = data_dir / "raw" / "alphafold" / "AF-P04637-F1-model_v6.pdb"

    # Load recommendations
    recommendations = json.loads(
        (reports_dir / "tables" / "tiered_recommendations.json").read_text()
    )

    print("=" * 80)
    print("DNA Binding Close-Up Images")
    print("=" * 80)
    print("\nGenerating 3 images per target:")
    print("  1. Wild-type (normal binding)")
    print("  2. Cancer mutant (disrupted binding)")
    print("  3. Rescued (restored binding)")
    print()

    for target in ["R175H", "R248Q", "R273H"]:
        if target not in recommendations:
            continue

        print(f"{target}:")

        modeled_dir = data_dir / "processed" / "modeled_structures" / target

        # Find structures
        cancer_pdb = modeled_dir / "cancer_build" / f"{target}_cancer.pdb"
        rescue_pdb = modeled_dir / "rescue_build" / f"{target}_rescued.pdb"

        if not cancer_pdb.exists() or not rescue_pdb.exists():
            print(f"  ✗ Modeled structures not found")
            continue

        # Image 1: Wild-type
        wt_img = closeups_dir / f"{target}_wildtype_binding.png"
        wt_script = scripts_dir / f"closeup_{target}_wildtype.pml"
        generate_closeup_script(
            wt_pdb, dna_pdb, target, "wildtype", wt_img, wt_script
        )
        print(f"  ✓ Wild-type: closeup_{target}_wildtype.pml")

        # Image 2: Cancer mutant
        cancer_img = closeups_dir / f"{target}_cancer_binding.png"
        cancer_script = scripts_dir / f"closeup_{target}_cancer.pml"
        generate_closeup_script(
            cancer_pdb, dna_pdb, target, "cancer", cancer_img, cancer_script
        )
        print(f"  ✓ Cancer: closeup_{target}_cancer.pml")

        # Image 3: Rescued
        rescue_img = closeups_dir / f"{target}_rescued_binding.png"
        rescue_script = scripts_dir / f"closeup_{target}_rescued.pml"
        generate_closeup_script(
            rescue_pdb, dna_pdb, target, "rescued", rescue_img, rescue_script
        )
        print(f"  ✓ Rescued: closeup_{target}_rescued.pml")

    print("\n" + "=" * 80)
    print("Close-up scripts generated!")
    print("=" * 80)
    print("\nRender all:")
    print("  for script in reports/pymol_scripts/closeup_*.pml; do")
    print("      /Applications/PyMOL.app/Contents/MacOS/PyMOL -c \"$script\"")
    print("  done")
    print(f"\nImages will be saved to: {closeups_dir}")


if __name__ == "__main__":
    sys.exit(main())
