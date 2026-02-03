#!/usr/bin/env python3
"""Create side-by-side structural comparison images showing REAL conformational changes.

Uses actual modeled mutant structures (from model_rescue_structures.py) to show:
- Wild-type: stable structure
- Cancer mutant: shows destabilization (actual structural distortion)
- Rescued: shows re-stabilization (structure restored)
"""
from __future__ import annotations

import sys
import json
from pathlib import Path

repo_root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(repo_root))


def generate_comparison_script(
    wt_pdb: Path,
    cancer_pdb: Path,
    rescue_pdb: Path,
    target: str,
    rescue_muts: list[str],
    output_png: Path,
    script_path: Path,
) -> None:
    """Generate PyMOL script for high-quality side-by-side comparison.

    Shows REAL structural differences between the three states.
    """
    import re

    # Parse cancer mutation
    match = re.match(r'([A-Z*])(\d+)([A-Z*])', target)
    if not match:
        raise ValueError(f"Invalid target: {target}")

    ref_aa, target_pos, mut_aa = match.groups()
    target_pos = int(target_pos)

    # Parse rescue positions
    rescue_positions = []
    for mut in rescue_muts:
        m = re.match(r'([A-Z*])(\d+)([A-Z*])', mut)
        if m:
            rescue_positions.append(int(m.group(2)))

    commands = [
        "# PyMOL Structural Comparison: Real Conformational Changes",
        f"# Target: {target}, Rescue: {', '.join(rescue_muts)}",
        "",
        "# Settings for high-quality render",
        "set ray_trace_mode, 1",
        "set ray_shadows, 1",
        "set antialias, 2",
        "set ambient, 0.5",
        "set specular, 0.6",
        "set depth_cue, 1",
        "set bg_rgb, [1, 1, 1]",
        "set cartoon_fancy_helices, 1",
        "set cartoon_smooth_loops, 0.8",
        "",
        "# Load all three structures",
        f"load {wt_pdb}, wt",
        f"load {cancer_pdb}, cancer",
        f"load {rescue_pdb}, rescued",
        "",
        "# Align all to wild-type for comparison",
        "align cancer, wt",
        "align rescued, wt",
        "",
        "# Position structures side-by-side",
        "translate [-50, 0, 0], wt",
        "# cancer stays at origin",
        "translate [50, 0, 0], rescued",
        "",
        "# Hide everything initially",
        "hide everything",
        "",
        "# === WILD-TYPE (Left) ===",
        "show cartoon, wt",
        "color gray70, wt",
        "",
        "# Highlight the position that WILL be mutated (for comparison)",
        f"select wt_future_mut, wt and resi {target_pos}",
        "show sticks, wt_future_mut",
        "color cyan, wt_future_mut",
        "set stick_radius, 0.25, wt_future_mut",
        "",
        "# Show stable structure indicator (green surface)",
        f"select wt_region, wt and (resi {target_pos} around 8)",
        "show surface, wt_region",
        "set transparency, 0.7, wt_region",
        "color palegreen, wt_region",
        "",
        "# === CANCER MUTANT (Middle) ===",
        "show cartoon, cancer",
        "color gray70, cancer",
        "",
        "# Highlight cancer mutation site (RED)",
        f"select cancer_mut, cancer and resi {target_pos}",
        "show spheres, cancer_mut",
        "color red, cancer_mut",
        "set sphere_scale, 1.8, cancer_mut",
        f"label cancer_mut and name CA, '{target} ({ref_aa}→{mut_aa})'",
        "",
        "# Show destabilized region (red surface)",
        f"select cancer_destab, cancer and (resi {target_pos} around 8)",
        "show surface, cancer_destab",
        "set transparency, 0.6, cancer_destab",
        "color salmon, cancer_destab",
        "",
        "# Show any structural distortions as sticks",
        "select cancer_distorted, cancer and (resi " + "+".join([str(p) for p in [target_pos] + list(range(target_pos-5, target_pos+6))]) + ")",
        "show sticks, cancer_distorted",
        "color orange, cancer_distorted",
        "set stick_radius, 0.2, cancer_distorted",
        "set stick_transparency, 0.5, cancer_distorted",
        "",
        "# === RESCUED (Right) ===",
        "show cartoon, rescued",
        "color gray70, rescued",
        "",
        "# Highlight rescue mutations (GREEN)",
    ]

    for i, (mut, pos) in enumerate(zip(rescue_muts, rescue_positions), 1):
        commands.extend([
            f"select rescue_mut_{i}, rescued and resi {pos}",
            f"show spheres, rescue_mut_{i}",
            f"color green, rescue_mut_{i}",
            f"set sphere_scale, 1.5, rescue_mut_{i}",
            f"label rescue_mut_{i} and name CA, '{mut}'",
        ])

    # Also show the cancer mutation site in rescued structure
    commands.extend([
        "",
        "# Show cancer mutation still present (but stabilized)",
        f"select rescued_cancer, rescued and resi {target_pos}",
        "show spheres, rescued_cancer",
        "color yellow, rescued_cancer",
        "set sphere_scale, 1.4, rescued_cancer",
        f"label rescued_cancer and name CA, '{target}'",
        "",
        "# Show re-stabilized region (green surface)",
        f"select rescued_stable, rescued and (resi {target_pos} around 8)",
        "show surface, rescued_stable",
        "set transparency, 0.7, rescued_stable",
        "color lightgreen, rescued_stable",
        "",
        "# === CALCULATE RMSD (shows actual structural differences) ===",
        "align cancer, wt, cycles=0",
        "print('='*60)",
        "print('STRUCTURAL DIFFERENCES (RMSD):')",
        "print('='*60)",
        "",
        "# Reset positions for RMSD calc",
        "translate [50, 0, 0], wt",
        "translate [0, 0, 0], cancer",
        "translate [-50, 0, 0], rescued",
        "",
        "# RMSD: wild-type vs cancer (shows destabilization)",
        "align cancer, wt",
        "print('Wild-type → Cancer mutant RMSD (destabilization):')",
        "",
        "# RMSD: cancer vs rescued (shows rescue effect)",
        "align rescued, cancer",
        "print('Cancer mutant → Rescued RMSD (re-stabilization):')",
        "",
        "# RMSD: wild-type vs rescued (how close to original)",
        "align rescued, wt",
        "print('Wild-type → Rescued RMSD (recovery):')",
        "print('='*60)",
        "",
        "# Reposition for visualization",
        "translate [-50, 0, 0], wt",
        "translate [0, 0, 0], cancer",
        "translate [50, 0, 0], rescued",
        "",
        "# === LABELS ===",
        "set label_size, 35",
        "set label_color, black",
        "set label_bg_color, white",
        "set label_bg_transparency, 0.3",
        "",
        "# Title labels",
        "pseudoatom label_wt, pos=[-50, 30, 0]",
        "label label_wt, 'WILD-TYPE\\nStable'",
        "",
        "pseudoatom label_cancer, pos=[0, 30, 0]",
        f"label label_cancer, 'CANCER\\n{target} Destabilized'",
        "",
        "pseudoatom label_rescued, pos=[50, 30, 0]",
        "label label_rescued, 'RESCUED\\n" + "+".join(rescue_muts) + "'",
        "",
        "# === ORIENT AND RENDER ===",
        "orient all",
        "zoom all, 6",
        "",
        "# Set viewport for widescreen",
        "viewport 3600, 1200",
        "",
        "# Ray trace at high quality",
        "print('Rendering high-quality image (this takes ~2-3 min)...')",
        "ray 3600, 1200",
        f"png {output_png}, dpi=300",
        "",
        "print('='*60)",
        f"print('Image saved: {output_png}')",
        "print('='*60)",
    ])

    script_path.write_text("\n".join(commands))
    print(f"  ✓ Generated comparison script: {script_path.name}")


def main():
    """Generate comparison scripts for modeled structures."""

    print("=" * 80)
    print("Structural Comparison Image Generator")
    print("=" * 80)
    print()

    data_dir = repo_root / "Data"
    reports_dir = repo_root / "reports"
    structures_dir = data_dir / "processed" / "modeled_structures"
    figures_dir = reports_dir / "figures"
    comparisons_dir = figures_dir / "structural_comparisons"
    scripts_dir = reports_dir / "pymol_scripts"

    comparisons_dir.mkdir(parents=True, exist_ok=True)
    scripts_dir.mkdir(parents=True, exist_ok=True)

    # Load recommendations
    recommendations_path = reports_dir / "tables" / "tiered_recommendations.json"
    if not recommendations_path.exists():
        print(f"ERROR: Recommendations not found: {recommendations_path}")
        return 1

    recommendations = json.loads(recommendations_path.read_text())

    print("Generating comparison scripts:\n")

    generated = []

    for target in ["R175H", "R248Q", "R273H"]:
        if target not in recommendations:
            continue

        best_data = (
            recommendations[target].get("best_1mut") or
            recommendations[target].get("best_2mut") or
            recommendations[target].get("best_3mut")
        )
        if not best_data:
            continue

        rescue_muts = best_data["rescue_mutations"]

        print(f"{target} + {', '.join(rescue_muts)}:")

        # Check for modeled structures
        target_dir = structures_dir / target
        wt_pdb = target_dir / "wild_type.pdb"

        # Find cancer mutant (prefer relaxed version)
        cancer_pdb = target_dir / f"{target}_cancer_relaxed.pdb"
        if not cancer_pdb.exists():
            cancer_pdb = target_dir / f"{target}_cancer.pdb"
            if not cancer_pdb.exists():
                cancer_pdb = target_dir / "cancer_build" / f"{target}_cancer.pdb"

        # Find rescued structure (prefer relaxed)
        rescue_pdb = target_dir / f"{target}_rescued_relaxed.pdb"
        if not rescue_pdb.exists():
            rescue_pdb = target_dir / f"{target}_rescued.pdb"
            if not rescue_pdb.exists():
                rescue_pdb = target_dir / "rescue_build" / f"{target}_rescued.pdb"

        # Check all structures exist
        if not all([wt_pdb.exists(), cancer_pdb.exists(), rescue_pdb.exists()]):
            print(f"  ✗ Missing modeled structures for {target}")
            print(f"     Run: python scripts/model_rescue_structures.py")
            continue

        # Generate comparison script
        output_png = comparisons_dir / f"{target}_structural_comparison.png"
        script_path = scripts_dir / f"compare_structures_{target}.pml"

        generate_comparison_script(
            wt_pdb,
            cancer_pdb,
            rescue_pdb,
            target,
            rescue_muts,
            output_png,
            script_path,
        )

        generated.append(script_path)

    print(f"\n{'=' * 80}")
    print(f"Generated {len(generated)} comparison scripts")
    print("=" * 80)
    print()

    if generated:
        print("To render comparisons, run:")
        for script in generated:
            print(f"  pymol -c {script}")
        print()
        print("Or render all at once:")
        print(f"  for script in {scripts_dir}/compare_structures_*.pml; do")
        print(f"      pymol -c \"$script\"")
        print(f"  done")
        print()
        print(f"Output images will be in: {comparisons_dir}")
    else:
        print("No scripts generated - run model_rescue_structures.py first")

    return 0


if __name__ == "__main__":
    sys.exit(main())
