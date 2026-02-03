#!/usr/bin/env python3
"""Create simple before/after videos showing DNA binding interactions.

For each target:
- Video 1: Cancer mutant with disrupted DNA binding
- Video 2: Rescued mutant with restored DNA binding
"""
from __future__ import annotations

import sys
from pathlib import Path

repo_root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(repo_root))


def generate_interaction_video_script(
    protein_pdb: Path,
    dna_pdb: Path,
    target: str,
    is_cancer: bool,
    output_path: Path,
    script_path: Path,
) -> None:
    """Generate PyMOL script showing protein-DNA interaction (60 frames = 2 sec)."""

    import re
    match = re.match(r'([A-Z*])(\d+)([A-Z*])', target)
    if not match:
        raise ValueError(f"Invalid target: {target}")

    _, target_pos, _ = match.groups()
    target_pos = int(target_pos)

    # DNA binding interface residues
    interface_residues = "241+248+273+280"

    color = "red" if is_cancer else "green"
    label = "CANCER MUTANT (Disrupted Binding)" if is_cancer else "RESCUED (Restored Binding)"

    commands = [
        "# p53-DNA Interaction Video",
        f"# {'Cancer' if is_cancer else 'Rescued'} mutant: {target}",
        "",
        "# Settings",
        "viewport 1280, 720",
        "set bg_rgb, [1, 1, 1]",
        "set antialias, 1",
        "set cartoon_fancy_helices, 1",
        "",
        "# Load structures",
        f"load {protein_pdb}, protein",
        f"load {dna_pdb}, dna",
        "",
        "# Align DNA to protein",
        "align dna, protein",
        "",
        "# Style protein",
        "hide everything",
        "show cartoon, protein",
        "color gray70, protein",
        "",
        "# Show DNA as sticks",
        "show cartoon, dna",
        "color blue, dna",
        "",
        "# Highlight mutation site",
        f"select mut_site, protein and resi {target_pos}",
        "show spheres, mut_site",
        f"color {color}, mut_site",
        "set sphere_scale, 1.6, mut_site",
        "",
        "# Highlight DNA binding interface",
        f"select interface, protein and resi {interface_residues}",
        "show sticks, interface",
        "color cyan, interface",
        "set stick_radius, 0.3, interface",
        "",
        "# Show hydrogen bonds (key difference!)",
        "distance hbonds, interface, dna, 3.5",
        "hide labels, hbonds",
        f"color {color}, hbonds",
        "set dash_width, 3, hbonds" if not is_cancer else "set dash_width, 2, hbonds",
        "set dash_gap, 0.3, hbonds" if is_cancer else "set dash_gap, 0.1, hbonds",
        "",
        "# Add surface to show binding pocket",
        f"select pocket, protein and (interface around 8)",
        "show surface, pocket",
        "set transparency, 0.75, pocket",
        f"color {'salmon' if is_cancer else 'palegreen'}, pocket",
        "",
        "# Label",
        "pseudoatom label, pos=[0,40,0]",
        f"label label, '{label}'",
        "set label_size, 50",
        f"set label_color, {color}",
        "",
        "# Orient to show binding interface",
        "orient interface",
        "zoom interface, 6",
        "",
        "# Create 60-frame movie (2 seconds @ 30fps)",
        "mset 1 x60",
        "",
        "# Smooth 360-degree rotation",
        "mview store, 1",
        "turn y, 360",
        "mview store, 60",
        "mview interpolate",
        "",
        "# Render frames (no ray-tracing for speed)",
        "set cache_frames, 0",
        f"mpng {output_path.parent / 'frames' / output_path.stem}_",
        "",
        "print('='*60)",
        "print('Frames rendered!')",
        "print('Create video:')",
        f"print('cd {output_path.parent}')",
        f"print('ffmpeg -framerate 30 -i frames/{output_path.stem}_%04d.png -c:v libx264 -pix_fmt yuv420p {output_path.name}')",
        "print('='*60)",
    ]

    script_path.write_text("\n".join(commands))


def main():
    """Generate interaction video scripts."""

    import json

    data_dir = repo_root / "Data"
    reports_dir = repo_root / "reports"
    animations_dir = reports_dir / "figures" / "animations"
    scripts_dir = reports_dir / "pymol_scripts"

    animations_dir.mkdir(parents=True, exist_ok=True)
    (animations_dir / "frames").mkdir(exist_ok=True)

    # DNA structure
    dna_pdb = data_dir / "raw" / "experimental_pdbs" / "2OCJ_full.pdb"

    # Load recommendations
    recommendations = json.loads(
        (reports_dir / "tables" / "tiered_recommendations.json").read_text()
    )

    print("=" * 80)
    print("Protein-DNA Interaction Videos: Before/After Rescue")
    print("=" * 80)
    print("\nGenerating 2 videos per target:")
    print("  1. Cancer mutant (disrupted binding)")
    print("  2. Rescued mutant (restored binding)")
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

        # Video 1: Cancer mutant (bad interaction)
        cancer_mp4 = animations_dir / f"{target}_cancer_binding.mp4"
        cancer_script = scripts_dir / f"interact_{target}_cancer.pml"

        generate_interaction_video_script(
            cancer_pdb, dna_pdb, target, is_cancer=True,
            output_path=cancer_mp4, script_path=cancer_script
        )
        print(f"  ✓ Cancer binding: interact_{target}_cancer.pml")

        # Video 2: Rescued mutant (good interaction)
        rescue_mp4 = animations_dir / f"{target}_rescued_binding.mp4"
        rescue_script = scripts_dir / f"interact_{target}_rescued.pml"

        generate_interaction_video_script(
            rescue_pdb, dna_pdb, target, is_cancer=False,
            output_path=rescue_mp4, script_path=rescue_script
        )
        print(f"  ✓ Rescued binding: interact_{target}_rescued.pml")

    print("\n" + "=" * 80)
    print("Interaction video scripts generated!")
    print("=" * 80)
    print("\nRender all:")
    print("  for script in reports/pymol_scripts/interact_*.pml; do")
    print("      /Applications/PyMOL.app/Contents/MacOS/PyMOL -c \"$script\"")
    print("  done")
    print("\nThen run ffmpeg commands (shown after each render)")


if __name__ == "__main__":
    sys.exit(main())
