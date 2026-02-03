#!/usr/bin/env python3
"""Create overlay images showing all 3 structures aligned to highlight differences."""
from pathlib import Path

repo_root = Path(__file__).resolve().parents[1]


def generate_overlay_script(
    wt_pdb: Path,
    cancer_pdb: Path,
    rescue_pdb: Path,
    dna_pdb: Path,
    target: str,
    output_path: Path,
    script_path: Path,
) -> None:
    """Generate PyMOL script showing all 3 structures overlaid."""

    import re
    match = re.match(r'([A-Z*])(\d+)([A-Z*])', target)
    if not match:
        raise ValueError(f"Invalid target: {target}")

    _, target_pos, _ = match.groups()
    target_pos = int(target_pos)

    interface_residues = f"{target_pos}+241+248+273+280"

    commands = [
        f"# Overlay Comparison: {target} - All 3 States Aligned",
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
        "# Load all structures",
        f"load {wt_pdb}, wt",
        f"load {cancer_pdb}, cancer",
        f"load {rescue_pdb}, rescued",
        f"load {dna_pdb}, dna",
        "",
        "# Align all to wild-type",
        "align cancer, wt",
        "align rescued, wt",
        "align dna, wt",
        "",
        "# Hide everything initially",
        "hide everything",
        "",
        "# === Show DNA (reference) ===",
        "show cartoon, dna",
        "color gray80, dna",
        "set cartoon_transparency, 0.7, dna",
        "",
        "# === Show proteins as thin cartoons ===",
        "show cartoon, wt",
        "show cartoon, cancer",
        "show cartoon, rescued",
        "set cartoon_transparency, 0.8, wt",
        "set cartoon_transparency, 0.8, cancer",
        "set cartoon_transparency, 0.8, rescued",
        "color cyan, wt",
        "color red, cancer",
        "color green, rescued",
        "",
        "# === Highlight mutation sites as spheres ===",
        f"select wt_site, wt and resi {target_pos}",
        f"select cancer_site, cancer and resi {target_pos}",
        f"select rescued_site, rescued and resi {target_pos}",
        "",
        "show spheres, wt_site",
        "show spheres, cancer_site",
        "show spheres, rescued_site",
        "",
        "color cyan, wt_site",
        "color red, cancer_site",
        "color yellow, rescued_site",
        "",
        "set sphere_scale, 1.5, wt_site",
        "set sphere_scale, 1.5, cancer_site",
        "set sphere_scale, 1.5, rescued_site",
        "",
        "# === Show binding interface ===",
        f"select wt_interface, wt and resi {interface_residues}",
        f"select cancer_interface, cancer and resi {interface_residues}",
        f"select rescued_interface, rescued and resi {interface_residues}",
        "",
        "show sticks, wt_interface",
        "show sticks, cancer_interface",
        "show sticks, rescued_interface",
        "",
        "set stick_radius, 0.25, wt_interface",
        "set stick_radius, 0.25, cancer_interface",
        "set stick_radius, 0.25, rescued_interface",
        "",
        "set stick_transparency, 0.3, wt_interface",
        "set stick_transparency, 0.3, cancer_interface",
        "set stick_transparency, 0.3, rescued_interface",
        "",
        "color cyan, wt_interface and elem C",
        "color red, cancer_interface and elem C",
        "color green, rescued_interface and elem C",
        "",
        "# === Show structural differences ===",
        "# Calculate RMSD to show differences",
        "print('='*60)",
        "print('STRUCTURAL ALIGNMENT:')",
        "align cancer, wt",
        "print('Cancer vs Wild-type RMSD above')",
        "align rescued, wt",
        "print('Rescued vs Wild-type RMSD above')",
        "print('='*60)",
        "",
        "# === Orient on binding interface ===",
        "center wt_interface",
        "orient wt_interface",
        "zoom wt_interface, 5",
        "turn y, 15",
        "turn x, 10",
        "",
        "# === Render ===",
        "viewport 1920, 1080",
        "ray 1920, 1080",
        f"png {output_path}, dpi=300",
        "",
        f"print('Overlay image saved: {output_path}')",
    ]

    script_path.write_text("\n".join(commands))
    print(f"  ✓ Overlay script: {script_path.name}")


def main():
    import json

    data_dir = repo_root / "Data"
    reports_dir = repo_root / "reports"
    figures_dir = reports_dir / "figures"
    overlays_dir = figures_dir / "overlay_comparisons"
    scripts_dir = reports_dir / "pymol_scripts"

    overlays_dir.mkdir(parents=True, exist_ok=True)

    # DNA structure
    dna_pdb = data_dir / "raw" / "experimental_pdbs" / "2OCJ_full.pdb"
    wt_pdb = data_dir / "raw" / "alphafold" / "AF-P04637-F1-model_v6.pdb"

    # Load recommendations
    recommendations = json.loads(
        (reports_dir / "tables" / "tiered_recommendations.json").read_text()
    )

    print("=" * 80)
    print("Creating Overlay Comparison Images")
    print("=" * 80)
    print("\nAll 3 structures aligned showing structural differences")
    print()

    for target in ["R175H", "R248Q", "R273H"]:
        if target not in recommendations:
            continue

        print(f"{target}:")

        modeled_dir = data_dir / "processed" / "modeled_structures" / target

        cancer_pdb = modeled_dir / "cancer_build" / f"{target}_cancer.pdb"
        rescue_pdb = modeled_dir / "rescue_build" / f"{target}_rescued.pdb"

        if not cancer_pdb.exists() or not rescue_pdb.exists():
            print(f"  ✗ Structures not found")
            continue

        overlay_img = overlays_dir / f"{target}_overlay_comparison.png"
        overlay_script = scripts_dir / f"overlay_{target}.pml"

        generate_overlay_script(
            wt_pdb, cancer_pdb, rescue_pdb, dna_pdb,
            target, overlay_img, overlay_script
        )

    print("\n" + "=" * 80)
    print("Overlay scripts generated!")
    print("=" * 80)
    print("\nRender:")
    print("  for script in reports/pymol_scripts/overlay_*.pml; do")
    print("      /Applications/PyMOL.app/Contents/MacOS/PyMOL -c \"$script\"")
    print("  done")


if __name__ == "__main__":
    import sys
    sys.exit(main())
