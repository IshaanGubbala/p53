#!/usr/bin/env python3
"""Create before/after rescue mutation animations showing DNA binding interactions.

Generates compelling molecular animations for science fair showing:
- Wild-type p53 binding DNA
- Cancer mutant destabilized/failing to bind
- Rescue mutant restoring DNA binding
"""
from __future__ import annotations

import sys
import json
from pathlib import Path
from typing import Any

# Add repo root to path
repo_root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(repo_root))


def generate_before_after_animation_script(
    wt_pdb: Path,
    mutant_pdb: Path,
    rescue_pdb: Path,
    dna_pdb: Path,
    target: str,
    rescue_muts: list[str],
    output_path: Path,
    script_path: Path,
) -> None:
    """Generate PyMOL script for before/after animation with DNA binding.

    Args:
        wt_pdb: Wild-type p53 structure
        mutant_pdb: Cancer mutant structure
        rescue_pdb: Rescued structure
        dna_pdb: DNA structure (from 2OCJ or similar)
        target: Cancer mutation (e.g., R175H)
        rescue_muts: List of rescue mutations
        output_path: Output movie file path
        script_path: PyMOL script output path
    """
    import re

    # Parse mutation positions
    target_match = re.match(r'([A-Z*])(\d+)([A-Z*])', target)
    if not target_match:
        raise ValueError(f"Invalid target mutation: {target}")

    _, target_pos, _ = target_match.groups()
    target_pos = int(target_pos)

    rescue_positions = []
    for mut in rescue_muts:
        match = re.match(r'([A-Z*])(\d+)([A-Z*])', mut)
        if match:
            rescue_positions.append(int(match.group(2)))

    commands = [
        "# PyMOL Animation: Before/After Rescue with DNA Binding",
        f"# Target: {target}",
        f"# Rescue: {', '.join(rescue_muts)}",
        "",
        "# Set movie parameters",
        "viewport 1920, 1080",
        "set ray_trace_frames, 1",
        "set ray_trace_mode, 1",
        "set antialias, 2",
        "set ambient, 0.5",
        "set specular, 0.6",
        "set ray_shadows, 1",
        "set depth_cue, 1",
        "set bg_rgb, [1, 1, 1]",
        "",
        "# Load structures",
        f"load {wt_pdb}, wildtype",
        f"load {mutant_pdb}, cancer_mutant",
        f"load {rescue_pdb}, rescued",
        f"load {dna_pdb}, dna",
        "",
        "# Align all structures to wild-type",
        "align cancer_mutant, wildtype",
        "align rescued, wildtype",
        "align dna, wildtype",
        "",
        "# Initially hide all",
        "hide everything",
        "",
        "# =====================",
        "# SCENE 1: Wild-Type (Frames 1-60)",
        "# =====================",
        "mset 1 x180",  # 180 frames total (6 seconds at 30 fps)
        "frame 1",
        "",
        "# Show wild-type only",
        "show cartoon, wildtype",
        "show cartoon, dna",
        "color gray70, wildtype",
        "color blue, dna",
        "set cartoon_fancy_helices, 1",
        "set cartoon_smooth_loops, 0.5",
        "",
        "# DNA-binding interface",
        f"select dna_contact, wildtype and (resi {target_pos} or resi 241 or resi 248 or resi 273 or resi 280)",
        "show sticks, dna_contact",
        "color cyan, dna_contact",
        "set stick_radius, 0.2, dna_contact",
        "",
        "# Show hydrogen bonds between DNA and protein",
        "distance hbonds_wt, (dna_contact), (dna), 3.5",
        "hide labels, hbonds_wt",
        "color yellow, hbonds_wt",
        "set dash_width, 3, hbonds_wt",
        "set dash_gap, 0.2, hbonds_wt",
        "",
        "# Add text label",
        "set label_size, 40",
        "set label_color, black",
        "set label_bg_color, white",
        "set label_bg_transparency, 0.3",
        "pseudoatom label_wt, pos=[50,50,50]",
        "label label_wt, 'WILD-TYPE: Stable DNA Binding'",
        "",
        "# Orient view",
        "orient dna_contact",
        "zoom dna_contact, 8",
        "",
        "# Gentle rotation for Scene 1",
        "mview store, 1",
        "turn y, 15",
        "mview store, 30",
        "turn y, 15",
        "mview store, 60",
        "mview interpolate",
        "",
        "# =====================",
        "# SCENE 2: Cancer Mutant (Frames 61-120)",
        "# =====================",
        "frame 61",
        "",
        "# Transition: fade out wild-type, fade in cancer mutant",
        "hide cartoon, wildtype",
        "hide sticks, dna_contact",
        "delete hbonds_wt",
        "delete label_wt",
        "",
        "# Show cancer mutant",
        "show cartoon, cancer_mutant",
        "color gray70, cancer_mutant",
        "",
        "# Highlight mutation site in RED",
        f"select cancer_site, cancer_mutant and resi {target_pos}",
        "show spheres, cancer_site",
        "color red, cancer_site",
        "set sphere_scale, 1.8, cancer_site",
        "",
        "# Show disrupted DNA binding interface",
        f"select mutant_interface, cancer_mutant and (resi {target_pos} or resi 241 or resi 248 or resi 273 or resi 280)",
        "show sticks, mutant_interface",
        "color orange, mutant_interface",
        "set stick_transparency, 0.4, mutant_interface",
        "",
        "# Show BROKEN hydrogen bonds (fewer/weaker)",
        "distance hbonds_mutant, (mutant_interface), (dna), 3.5",
        "hide labels, hbonds_mutant",
        "color red, hbonds_mutant",
        "set dash_width, 2, hbonds_mutant",
        "set dash_gap, 0.4, hbonds_mutant",
        "",
        "# Structural distortion visualization",
        f"select distorted_region, cancer_mutant and (cancer_site around 8)",
        "show surface, distorted_region",
        "set transparency, 0.7, distorted_region",
        "color red, distorted_region",
        "",
        "# Label",
        "pseudoatom label_mutant, pos=[50,50,50]",
        f"label label_mutant, 'CANCER MUTANT ({target}): Destabilized, Poor DNA Binding'",
        "",
        "# Re-orient",
        "orient cancer_site",
        "zoom cancer_site, 8",
        "",
        "# Wobble motion to show instability",
        "mview store, 61",
        "turn y, -10",
        "turn x, 5",
        "mview store, 75",
        "turn y, 20",
        "turn x, -10",
        "mview store, 90",
        "turn y, -10",
        "turn x, 5",
        "mview store, 105",
        "turn y, 5",
        "mview store, 120",
        "mview interpolate",
        "",
        "# =====================",
        "# SCENE 3: Rescued Structure (Frames 121-180)",
        "# =====================",
        "frame 121",
        "",
        "# Transition: fade out cancer mutant, fade in rescued",
        "hide cartoon, cancer_mutant",
        "hide spheres, cancer_site",
        "hide sticks, mutant_interface",
        "hide surface, distorted_region",
        "delete hbonds_mutant",
        "delete label_mutant",
        "",
        "# Show rescued structure",
        "show cartoon, rescued",
        "color gray70, rescued",
        "",
        "# Highlight rescue mutations in GREEN",
    ]

    for i, (mut, pos) in enumerate(zip(rescue_muts, rescue_positions), 1):
        commands.extend([
            f"select rescue_{i}, rescued and resi {pos}",
            f"show spheres, rescue_{i}",
            f"color green, rescue_{i}",
            f"set sphere_scale, 1.5, rescue_{i}",
            f"label rescue_{i} and name CA, '{mut}'",
        ])

    # Highlight original cancer mutation site (now present in rescued structure)
    commands.extend([
        f"select rescued_cancer_site, rescued and resi {target_pos}",
        "show spheres, rescued_cancer_site",
        "color salmon, rescued_cancer_site",
        "set sphere_scale, 1.6, rescued_cancer_site",
        f"label rescued_cancer_site and name CA, '{target}'",
        "",
        "# Show RESTORED DNA binding interface",
        f"select rescued_interface, rescued and (resi {target_pos} or resi 241 or resi 248 or resi 273 or resi 280)",
        "show sticks, rescued_interface",
        "color cyan, rescued_interface",
        "set stick_radius, 0.2, rescued_interface",
        "",
        "# Show RESTORED hydrogen bonds",
        "distance hbonds_rescued, (rescued_interface), (dna), 3.5",
        "hide labels, hbonds_rescued",
        "color lime, hbonds_rescued",
        "set dash_width, 3, hbonds_rescued",
        "set dash_gap, 0.2, hbonds_rescued",
        "",
        "# Highlight stabilized region (green surface)",
        f"select stabilized_region, rescued and (rescued_cancer_site around 8)",
        "show surface, stabilized_region",
        "set transparency, 0.7, stabilized_region",
        "color green, stabilized_region",
        "",
        "# Label",
        "pseudoatom label_rescued, pos=[50,50,50]",
        f"label label_rescued, 'RESCUED: Stability Restored, DNA Binding Recovered'",
        "",
        "# Orient view",
        "orient rescued_interface",
        "zoom rescued_interface, 8",
        "",
        "# Smooth rotation for Scene 3",
        "mview store, 121",
        "turn y, 30",
        "mview store, 150",
        "turn y, 30",
        "mview store, 180",
        "mview interpolate",
        "",
        "# =====================",
        "# RENDER ANIMATION",
        "# =====================",
        "",
        "# Set output format and render settings",
        f"set movie_fps, 30",
        "set ray_trace_frames, 1",
        "set cache_frames, 0",
        "",
        "print('Animation rendering started (this takes ~10 min)...')",
        "# Use mpng to render individual PNG frames, then ffmpeg will combine",
        f"mpng {str(output_path.parent / 'frames' / output_path.stem)}_",
        "",
        "print('Frames rendered! Use ffmpeg to create video:')",
        f"print('cd {output_path.parent}')",
        f"print('ffmpeg -framerate 30 -i frames/{output_path.stem}_%04d.png -c:v libx264 -pix_fmt yuv420p {output_path.name}')",
        "",
    ])

    # Write script
    script_path.write_text("\n".join(commands))
    print(f"✓ Generated animation script: {script_path}")


def generate_side_by_side_comparison_script(
    wt_pdb: Path,
    mutant_pdb: Path,
    rescue_pdb: Path,
    dna_pdb: Path,
    target: str,
    rescue_muts: list[str],
    output_path: Path,
    script_path: Path,
) -> None:
    """Generate side-by-side comparison showing all three states simultaneously."""

    import re

    target_match = re.match(r'([A-Z*])(\d+)([A-Z*])', target)
    if not target_match:
        raise ValueError(f"Invalid target mutation: {target}")

    _, target_pos, _ = target_match.groups()
    target_pos = int(target_pos)

    rescue_positions = []
    for mut in rescue_muts:
        match = re.match(r'([A-Z*])(\d+)([A-Z*])', mut)
        if match:
            rescue_positions.append(int(match.group(2)))

    commands = [
        "# PyMOL Side-by-Side Comparison",
        f"# Target: {target}, Rescue: {', '.join(rescue_muts)}",
        "",
        "# Load structures",
        f"load {wt_pdb}, wt",
        f"load {mutant_pdb}, mutant",
        f"load {rescue_pdb}, rescued",
        f"load {dna_pdb}, dna_wt",
        f"load {dna_pdb}, dna_mut",
        f"load {dna_pdb}, dna_res",
        "",
        "# Align all",
        "align mutant, wt",
        "align rescued, wt",
        "align dna_wt, wt",
        "align dna_mut, mutant",
        "align dna_res, rescued",
        "",
        "# Position structures side-by-side",
        "translate [-40, 0, 0], wt",
        "translate [-40, 0, 0], dna_wt",
        "translate [40, 0, 0], rescued",
        "translate [40, 0, 0], dna_res",
        "# mutant stays at origin",
        "",
        "# Style all proteins",
        "hide everything",
        "show cartoon, wt or mutant or rescued",
        "show cartoon, dna_*",
        "color gray70, wt or mutant or rescued",
        "color blue, dna_*",
        "set cartoon_fancy_helices, 1",
        "",
        "# Wild-type: Cyan highlights",
        f"select wt_interface, wt and resi {target_pos}+241+248+273+280",
        "show sticks, wt_interface",
        "color cyan, wt_interface",
        "distance hb_wt, wt_interface, dna_wt, 3.5",
        "hide labels, hb_wt",
        "color lime, hb_wt",
        "set dash_width, 3, hb_wt",
        "",
        "# Cancer mutant: Red highlights",
        f"select mut_site, mutant and resi {target_pos}",
        "show spheres, mut_site",
        "color red, mut_site",
        "set sphere_scale, 1.8, mut_site",
        "distance hb_mut, mutant, dna_mut, 3.5",
        "hide labels, hb_mut",
        "color red, hb_mut",
        "set dash_width, 2, hb_mut",
        "set dash_gap, 0.5, hb_mut",
        "",
        "# Rescued: Green highlights",
        f"select res_site, rescued and resi {target_pos}",
        "show spheres, res_site",
        "color salmon, res_site",
        "set sphere_scale, 1.6, res_site",
    ]

    for i, pos in enumerate(rescue_positions, 1):
        commands.extend([
            f"select rescue_{i}, rescued and resi {pos}",
            f"show spheres, rescue_{i}",
            f"color green, rescue_{i}",
            f"set sphere_scale, 1.4, rescue_{i}",
        ])

    commands.extend([
        "distance hb_res, rescued, dna_res, 3.5",
        "hide labels, hb_res",
        "color lime, hb_res",
        "set dash_width, 3, hb_res",
        "",
        "# Labels",
        "pseudoatom label1, pos=[-40, 25, 0]",
        "label label1, 'WILD-TYPE\\nStable DNA Binding'",
        "pseudoatom label2, pos=[0, 25, 0]",
        f"label label2, 'CANCER ({target})\\nDestabilized'",
        "pseudoatom label3, pos=[40, 25, 0]",
        "label label3, 'RESCUED\\n" + " + ".join(rescue_muts) + "'",
        "",
        "set label_size, 30",
        "set label_color, black",
        "set label_bg_color, white",
        "set label_bg_transparency, 0.3",
        "",
        "# Orient and render",
        "orient all",
        "zoom all, 5",
        "viewport 3600, 1200",
        "bg_color white",
        "set ray_shadows, 1",
        "set antialias, 2",
        "ray 3600, 1200",
        f"png {output_path}, dpi=300",
        "",
        "print('Side-by-side comparison complete!')",
    ])

    script_path.write_text("\n".join(commands))
    print(f"✓ Generated side-by-side script: {script_path}")


def main():
    """Generate before/after rescue animation scripts."""

    print("=" * 80)
    print("p53 Rescue Animation Generator: Before/After with DNA Binding")
    print("=" * 80)

    data_dir = repo_root / "Data"
    reports_dir = repo_root / "reports"
    figures_dir = reports_dir / "figures"
    animations_dir = figures_dir / "animations"
    scripts_dir = reports_dir / "pymol_scripts"

    animations_dir.mkdir(parents=True, exist_ok=True)
    scripts_dir.mkdir(parents=True, exist_ok=True)

    # Find structures
    alphafold_dir = data_dir / "raw" / "alphafold"
    experimental_dir = data_dir / "raw" / "experimental_pdbs"

    # Wild-type structure
    wt_pdb = alphafold_dir / "AF-P04637-F1-model_v6.pdb"
    if not wt_pdb.exists():
        pdb_files = list(alphafold_dir.glob("*.pdb"))
        if pdb_files:
            wt_pdb = pdb_files[0]
        else:
            print(f"ERROR: No wild-type PDB found in {alphafold_dir}")
            return 1

    # DNA structure (2OCJ full complex with DNA)
    dna_pdb = experimental_dir / "2OCJ_full.pdb"
    if not dna_pdb.exists():
        # Try chain-only version
        dna_pdb = experimental_dir / "2OCJ_chainA.pdb"
        if not dna_pdb.exists():
            print(f"WARNING: DNA structure not found")
            print("         Animation will show protein only (no DNA binding)")
            # Use wild-type as placeholder
            dna_pdb = wt_pdb

    print(f"\nStructures:")
    print(f"  Wild-type: {wt_pdb.name}")
    print(f"  DNA: {dna_pdb.name}")

    # Load rescue recommendations
    recommendations_path = reports_dir / "tables" / "tiered_recommendations.json"
    if not recommendations_path.exists():
        print(f"\nERROR: Recommendations not found at {recommendations_path}")
        print("       Run 'python -m experiments.run_make_report' first")
        return 1

    recommendations = json.loads(recommendations_path.read_text())

    print(f"\n{'=' * 80}")
    print("Generating animation scripts for top rescue candidates:")
    print("=" * 80)

    # Generate animations for best rescue of each target
    for target in ["R175H", "R248Q", "R273H"]:
        if target not in recommendations:
            continue

        best_data = recommendations[target].get("best_1mut") or recommendations[target].get("best_2mut") or recommendations[target].get("best_3mut")
        if not best_data:
            print(f"  ✗ No rescue data found for {target}")
            continue

        rescue_muts = best_data["rescue_mutations"]

        print(f"\n{target} + {', '.join(rescue_muts)}:")

        # Use REAL modeled structures from EvoEF2
        modeled_dir = data_dir / "processed" / "modeled_structures" / target

        # Find cancer mutant structure (prefer relaxed if MD was run)
        mutant_pdb = modeled_dir / f"{target}_cancer_relaxed.pdb"
        if not mutant_pdb.exists():
            mutant_pdb = modeled_dir / "cancer_build" / f"{target}_cancer.pdb"

        # Find rescued structure (prefer relaxed)
        rescue_pdb = modeled_dir / f"{target}_rescued_relaxed.pdb"
        if not rescue_pdb.exists():
            rescue_pdb = modeled_dir / "rescue_build" / f"{target}_rescued.pdb"

        # Check if modeled structures exist
        if not mutant_pdb.exists() or not rescue_pdb.exists():
            print(f"  ✗ Modeled structures not found for {target}")
            print(f"     Run: python scripts/model_rescue_structures.py")
            continue

        print(f"  Using REAL modeled structures:")
        print(f"    - Wild-type: {wt_pdb.name}")
        print(f"    - Cancer mutant: {mutant_pdb.name}")
        print(f"    - Rescued: {rescue_pdb.name}")

        # 1. Before/after sequential animation
        anim_output = animations_dir / f"rescue_animation_{target}.mp4"
        anim_script = scripts_dir / f"animate_rescue_{target}.pml"

        try:
            generate_before_after_animation_script(
                wt_pdb, mutant_pdb, rescue_pdb, dna_pdb,
                target, rescue_muts, anim_output, anim_script
            )
            print(f"  ✓ Animation script: {anim_script.name}")
        except Exception as e:
            print(f"  ✗ Animation script failed: {e}")

        # 2. Side-by-side comparison
        comp_output = animations_dir / f"rescue_comparison_{target}.png"
        comp_script = scripts_dir / f"compare_rescue_{target}.pml"

        try:
            generate_side_by_side_comparison_script(
                wt_pdb, mutant_pdb, rescue_pdb, dna_pdb,
                target, rescue_muts, comp_output, comp_script
            )
            print(f"  ✓ Comparison script: {comp_script.name}")
        except Exception as e:
            print(f"  ✗ Comparison script failed: {e}")

    print(f"\n{'=' * 80}")
    print("Scripts generated successfully!")
    print("=" * 80)

    print(f"\nNext steps:")
    print(f"  1. Execute PyMOL scripts:")
    print(f"     pymol -cq {scripts_dir}/animate_rescue_R175H.pml")
    print(f"  2. View animations:")
    print(f"     open {animations_dir}/")
    print(f"\n  NOTE: These scripts use wild-type PDB as placeholders.")
    print(f"        For real mutations, you need to:")
    print(f"        - Model cancer mutant structures (e.g., with FoldX or Rosetta)")
    print(f"        - Model rescue mutant structures")
    print(f"        - Use actual DNA-bound structure (e.g., 2OCJ_full.pdb with DNA)")

    return 0


if __name__ == "__main__":
    sys.exit(main())
