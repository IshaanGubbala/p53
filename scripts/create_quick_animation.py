#!/usr/bin/env python3
"""Create quick animations using real mutant structures (no ray-tracing for speed)."""
from __future__ import annotations

import sys
from pathlib import Path

repo_root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(repo_root))


def generate_quick_animation(
    wt_pdb: Path,
    cancer_pdb: Path,
    rescue_pdb: Path,
    target: str,
    output_path: Path,
    script_path: Path,
) -> None:
    """Generate fast PyMOL animation script (no ray-tracing)."""

    import re
    match = re.match(r'([A-Z*])(\d+)([A-Z*])', target)
    if not match:
        raise ValueError(f"Invalid target: {target}")

    _, target_pos, _ = match.groups()
    target_pos = int(target_pos)

    commands = [
        "# Quick p53 Rescue Animation (No Ray-Tracing for Speed)",
        f"# Target: {target}",
        "",
        "# Settings",
        "viewport 1280, 720",
        "set antialias, 1",
        "set bg_rgb, [1, 1, 1]",
        "set cartoon_fancy_helices, 1",
        "set cartoon_smooth_loops, 1",
        "",
        "# Load all three structures",
        f"load {wt_pdb}, wt",
        f"load {cancer_pdb}, cancer",
        f"load {rescue_pdb}, rescued",
        "",
        "# Align",
        "align cancer, wt",
        "align rescued, wt",
        "",
        "# Hide all initially",
        "hide everything",
        "",
        "# Create 90-frame movie (3 seconds @ 30fps)",
        "mset 1 x90",
        "",
        "# === SCENE 1: Wild-Type (Frames 1-30) ===",
        "frame 1",
        "show cartoon, wt",
        "color gray70, wt",
        f"select wt_site, wt and resi {target_pos}",
        "show spheres, wt_site",
        "color cyan, wt_site",
        "set sphere_scale, 1.5, wt_site",
        "",
        "# Create label",
        "pseudoatom label1, pos=[30,30,0]",
        "label label1, 'WILD-TYPE (Healthy)'",
        "set label_size, 40",
        "set label_color, black",
        "",
        "# Orient",
        "orient wt",
        "zoom wt, 5",
        "",
        "# Smooth rotation",
        "mview store, 1",
        "turn y, 60",
        "mview store, 30",
        "",
        "# === SCENE 2: Cancer Mutant (Frames 31-60) ===",
        "frame 31",
        "",
        "# Transition",
        "hide cartoon, wt",
        "hide spheres, wt_site",
        "delete label1",
        "",
        "# Show cancer mutant",
        "show cartoon, cancer",
        "color gray70, cancer",
        f"select cancer_site, cancer and resi {target_pos}",
        "show spheres, cancer_site",
        "color red, cancer_site",
        "set sphere_scale, 1.8, cancer_site",
        "",
        "# Show destabilized region",
        f"select destab, cancer and (cancer_site around 8)",
        "show surface, destab",
        "set transparency, 0.7, destab",
        "color salmon, destab",
        "",
        "# Label",
        "pseudoatom label2, pos=[30,30,0]",
        f"label label2, 'CANCER MUTANT ({target})'",
        "set label_color, red",
        "",
        "# Orient",
        "orient cancer",
        "zoom cancer, 5",
        "",
        "# Smooth rotation",
        "mview store, 31",
        "turn y, 60",
        "mview store, 60",
        "",
        "# === SCENE 3: Rescued (Frames 61-90) ===",
        "frame 61",
        "",
        "# Transition",
        "hide cartoon, cancer",
        "hide spheres, cancer_site",
        "hide surface, destab",
        "delete label2",
        "",
        "# Show rescued",
        "show cartoon, rescued",
        "color gray70, rescued",
        f"select rescued_site, rescued and resi {target_pos}",
        "show spheres, rescued_site",
        "color yellow, rescued_site",
        "set sphere_scale, 1.6, rescued_site",
        "",
        "# Show re-stabilized region",
        f"select stab, rescued and (rescued_site around 8)",
        "show surface, stab",
        "set transparency, 0.7, stab",
        "color palegreen, stab",
        "",
        "# Label",
        "pseudoatom label3, pos=[30,30,0]",
        "label label3, 'RESCUED (Stability Restored)'",
        "set label_color, green",
        "",
        "# Orient",
        "orient rescued",
        "zoom rescued, 5",
        "",
        "# Smooth rotation",
        "mview store, 61",
        "turn y, 60",
        "mview store, 90",
        "",
        "# Interpolate all views",
        "mview interpolate",
        "mview reinterpolate",
        "",
        "# === RENDER ===",
        f"# Render to PNG sequence (NO ray-tracing for speed)",
        "set cache_frames, 0",
        f"mpng {output_path.parent / 'frames' / output_path.stem}_",
        "",
        "print('='*60)",
        "print('Animation frames saved!')",
        "print('Create video with ffmpeg:')",
        f"print('cd {output_path.parent}')",
        f"print('ffmpeg -framerate 30 -i frames/{output_path.stem}_%04d.png -vf scale=1280:720 -c:v libx264 -pix_fmt yuv420p {output_path.name}')",
        "print('='*60)",
    ]

    script_path.write_text("\n".join(commands))
    print(f"  ✓ Quick animation script: {script_path.name}")


def main():
    """Generate quick animation scripts."""

    import json

    data_dir = repo_root / "Data"
    reports_dir = repo_root / "reports"
    animations_dir = reports_dir / "figures" / "animations"
    scripts_dir = reports_dir / "pymol_scripts"

    animations_dir.mkdir(parents=True, exist_ok=True)
    (animations_dir / "frames").mkdir(exist_ok=True)
    scripts_dir.mkdir(parents=True, exist_ok=True)

    # Wild-type structure
    wt_pdb = data_dir / "raw" / "alphafold" / "AF-P04637-F1-model_v6.pdb"

    # Load recommendations
    recommendations = json.loads(
        (reports_dir / "tables" / "tiered_recommendations.json").read_text()
    )

    print("=" * 80)
    print("Quick Animation Generator: Real Structural Changes")
    print("=" * 80)

    for target in ["R175H", "R248Q", "R273H"]:
        if target not in recommendations:
            continue

        print(f"\n{target}:")

        modeled_dir = data_dir / "processed" / "modeled_structures" / target

        # Find structures
        cancer_pdb = modeled_dir / "cancer_build" / f"{target}_cancer.pdb"
        rescue_pdb = modeled_dir / "rescue_build" / f"{target}_rescued.pdb"

        if not cancer_pdb.exists() or not rescue_pdb.exists():
            print(f"  ✗ Modeled structures not found")
            continue

        output_mp4 = animations_dir / f"quick_rescue_{target}.mp4"
        script_pml = scripts_dir / f"quick_animate_{target}.pml"

        generate_quick_animation(
            wt_pdb, cancer_pdb, rescue_pdb, target, output_mp4, script_pml
        )

    print("\n" + "=" * 80)
    print("Quick animation scripts generated!")
    print("=" * 80)
    print("\nTo render:")
    print("  for script in reports/pymol_scripts/quick_animate_*.pml; do")
    print("      /Applications/PyMOL.app/Contents/MacOS/PyMOL -c \"$script\"")
    print("  done")
    print("\nThen create videos with ffmpeg (see output after rendering)")


if __name__ == "__main__":
    sys.exit(main())
