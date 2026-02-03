#!/usr/bin/env python3
"""Model actual mutant structures showing real conformational changes.

Uses EvoEF2 BuildMutant to create:
1. Wild-type structure (baseline)
2. Cancer mutant structure (R175H) - shows destabilization
3. Rescued structure (R175H + M133L) - shows re-stabilization

This gives REAL structural differences for animation, not just color changes.
"""
from __future__ import annotations

import sys
import json
import shutil
import subprocess
from pathlib import Path
from typing import Any

# Add repo root to path
repo_root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(repo_root))


def parse_mutation(mut_str: str) -> tuple[str, int, str]:
    """Parse mutation string like 'R175H' into (ref, pos, alt)."""
    import re
    match = re.match(r'([A-Z*])(\d+)([A-Z*])', mut_str)
    if not match:
        raise ValueError(f"Invalid mutation: {mut_str}")
    ref, pos, alt = match.groups()
    return ref, int(pos), alt


def create_mutation_file(mutations: list[str], output_path: Path, chain_id: str = "A") -> None:
    """Create individual_list.txt for EvoEF2 BuildMutant.

    Format (with chain ID):
    RA175H;  (R at chain A position 175 -> H)
    MA133L;  (M at chain A position 133 -> L)
    """
    # Convert mutations like "R175H" to "RA175H" (add chain ID)
    formatted_muts = []
    for mut in mutations:
        ref, pos, alt = parse_mutation(mut)
        formatted_mut = f"{ref}{chain_id}{pos}{alt}"
        formatted_muts.append(formatted_mut)

    lines = [f"{mut};" for mut in formatted_muts]
    output_path.write_text("\n".join(lines) + "\n")
    print(f"  ✓ Created mutation file: {output_path}")


def run_evoef2_buildmutant(
    evoef2_binary: Path,
    wild_pdb: Path,
    mutation_file: Path,
    output_dir: Path,
    mutant_name: str,
) -> Path | None:
    """Run EvoEF2 BuildMutant to create mutant structure.

    Args:
        evoef2_binary: Path to EvoEF2 executable
        wild_pdb: Wild-type PDB file
        mutation_file: File with mutation list
        output_dir: Output directory
        mutant_name: Name for output (e.g., "R175H_cancer")

    Returns:
        Path to generated mutant PDB, or None if failed
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # EvoEF2 BuildMutant command
    cmd = [
        str(evoef2_binary),
        "--command=BuildMutant",
        f"--pdb={wild_pdb}",
        f"--mutant_file={mutation_file}",  # Note: underscore not hyphen
    ]

    print(f"  → Running EvoEF2 BuildMutant for {mutant_name}...")
    print(f"     Command: {' '.join(cmd)}")

    try:
        result = subprocess.run(
            cmd,
            cwd=output_dir,
            capture_output=True,
            text=True,
            timeout=300,  # 5 minutes max
        )

        if result.returncode != 0:
            print(f"  ✗ EvoEF2 failed:")
            print(f"     {result.stderr[:200]}")
            return None

        # EvoEF2 creates files like "*_Model_0001.pdb" in the output directory
        potential_files = list(output_dir.glob("*Model_0001.pdb"))

        if not potential_files:
            # Try alternate patterns
            potential_files = list(output_dir.glob("*Model_*.pdb"))

        if not potential_files:
            print(f"  ✗ No mutant PDB found in {output_dir}")
            return None

        mutant_pdb = sorted(potential_files)[0]

        # Rename to meaningful name
        final_pdb = output_dir / f"{mutant_name}.pdb"
        shutil.copy(mutant_pdb, final_pdb)

        print(f"  ✓ Created mutant structure: {final_pdb}")
        return final_pdb

    except subprocess.TimeoutExpired:
        print(f"  ✗ EvoEF2 timed out after 5 minutes")
        return None
    except Exception as e:
        print(f"  ✗ EvoEF2 error: {e}")
        return None


def run_molecular_dynamics_relaxation(
    mutant_pdb: Path,
    output_pdb: Path,
    steps: int = 1000,
) -> Path | None:
    """Run short MD relaxation to show conformational change.

    Uses OpenMM for quick energy minimization + short MD.
    This shows how the mutant structure actually relaxes/distorts.

    Args:
        mutant_pdb: Input mutant structure
        output_pdb: Output relaxed structure
        steps: Number of MD steps (1000 = ~1 picosecond)

    Returns:
        Path to relaxed structure, or None if failed
    """
    try:
        from openmm import app
        from openmm import unit
        import openmm as mm
    except ImportError:
        print("  ⚠ OpenMM not installed - skipping MD relaxation")
        print("     Install with: conda install -c conda-forge openmm")
        return None

    print(f"  → Running MD relaxation ({steps} steps)...")

    try:
        # Load structure
        pdb = app.PDBFile(str(mutant_pdb))

        # Set up force field
        forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

        # Create system
        system = forcefield.createSystem(
            pdb.topology,
            nonbondedMethod=app.NoCutoff,
            constraints=app.HBonds,
        )

        # Energy minimization first
        integrator = mm.LangevinMiddleIntegrator(
            300*unit.kelvin,
            1/unit.picosecond,
            0.002*unit.picoseconds,
        )
        simulation = app.Simulation(pdb.topology, system, integrator)
        simulation.context.setPositions(pdb.positions)

        print(f"     Minimizing energy...")
        simulation.minimizeEnergy(maxIterations=500)

        # Short MD to relax structure
        print(f"     Running {steps} MD steps...")
        simulation.step(steps)

        # Save final structure
        state = simulation.context.getState(getPositions=True)
        with open(output_pdb, 'w') as f:
            app.PDBFile.writeFile(
                simulation.topology,
                state.getPositions(),
                f,
            )

        print(f"  ✓ Relaxed structure saved: {output_pdb}")
        return output_pdb

    except Exception as e:
        print(f"  ✗ MD relaxation failed: {e}")
        return None


def main():
    """Generate actual mutant structures for all rescue targets."""

    print("=" * 80)
    print("p53 Mutant Structure Modeling: Real Conformational Changes")
    print("=" * 80)
    print()

    data_dir = repo_root / "Data"
    reports_dir = repo_root / "reports"
    structures_dir = data_dir / "processed" / "modeled_structures"
    structures_dir.mkdir(parents=True, exist_ok=True)

    # Find wild-type structure
    wt_pdb = data_dir / "raw" / "alphafold" / "AF-P04637-F1-model_v6.pdb"
    if not wt_pdb.exists():
        print(f"ERROR: Wild-type PDB not found: {wt_pdb}")
        return 1

    print(f"Wild-type structure: {wt_pdb.name}\n")

    # Find EvoEF2 binary
    evoef2_binary = shutil.which("EvoEF2")
    if not evoef2_binary:
        # Try common locations
        possible_paths = [
            Path("/Users/ishaangubbala/Downloads/EvoEF2-master/EvoEF2"),  # From config
            repo_root / "bin" / "EvoEF2",
            Path.home() / "bin" / "EvoEF2",
            Path("/usr/local/bin/EvoEF2"),
        ]
        for path in possible_paths:
            if path.exists():
                evoef2_binary = path
                break

    if not evoef2_binary:
        print("ERROR: EvoEF2 not found!")
        print("  Install from: https://github.com/tommyhuangthu/EvoEF2")
        print("  Or specify path in this script")
        return 1

    evoef2_binary = Path(evoef2_binary)
    print(f"EvoEF2 binary: {evoef2_binary}\n")

    # Load rescue recommendations
    recommendations_path = reports_dir / "tables" / "tiered_recommendations.json"
    if not recommendations_path.exists():
        print(f"ERROR: Recommendations not found: {recommendations_path}")
        print("  Run 'python -m experiments.run_make_report' first")
        return 1

    recommendations = json.loads(recommendations_path.read_text())

    print("=" * 80)
    print("Modeling Structures:")
    print("=" * 80)
    print()

    # Process each target
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

        print(f"{'=' * 80}")
        print(f"Target: {target} + {', '.join(rescue_muts)}")
        print(f"{'=' * 80}")

        target_dir = structures_dir / target
        target_dir.mkdir(exist_ok=True)

        # 1. Copy wild-type (baseline)
        wt_copy = target_dir / "wild_type.pdb"
        shutil.copy(wt_pdb, wt_copy)
        print(f"1. Wild-type: {wt_copy.name}")

        # 2. Model cancer mutant (e.g., R175H)
        print(f"\n2. Cancer mutant ({target}):")
        cancer_mut_file = target_dir / "cancer_mutations.txt"
        create_mutation_file([target], cancer_mut_file)

        cancer_pdb = run_evoef2_buildmutant(
            evoef2_binary,
            wt_pdb,
            cancer_mut_file,
            target_dir / "cancer_build",
            f"{target}_cancer"
        )

        # Optional: Relax cancer mutant with MD to show destabilization
        if cancer_pdb:
            cancer_relaxed = target_dir / f"{target}_cancer_relaxed.pdb"
            relaxed = run_molecular_dynamics_relaxation(
                cancer_pdb,
                cancer_relaxed,
                steps=1000,
            )
            if relaxed:
                cancer_pdb = relaxed

        # 3. Model rescued structure (cancer + rescue mutations)
        print(f"\n3. Rescued structure ({target} + {', '.join(rescue_muts)}):")
        all_mutations = [target] + rescue_muts
        rescue_mut_file = target_dir / "rescue_mutations.txt"
        create_mutation_file(all_mutations, rescue_mut_file)

        rescue_pdb = run_evoef2_buildmutant(
            evoef2_binary,
            wt_pdb,
            rescue_mut_file,
            target_dir / "rescue_build",
            f"{target}_rescued"
        )

        # Optional: Relax rescued structure with MD
        if rescue_pdb:
            rescue_relaxed = target_dir / f"{target}_rescued_relaxed.pdb"
            relaxed = run_molecular_dynamics_relaxation(
                rescue_pdb,
                rescue_relaxed,
                steps=1000,
            )
            if relaxed:
                rescue_pdb = relaxed

        print(f"\n{'=' * 80}")
        print(f"Summary for {target}:")
        print(f"  Wild-type:  {wt_copy}")
        if cancer_pdb:
            print(f"  Cancer:     {cancer_pdb}")
        else:
            print(f"  Cancer:     FAILED")
        if rescue_pdb:
            print(f"  Rescued:    {rescue_pdb}")
        else:
            print(f"  Rescued:    FAILED")
        print()

    print("=" * 80)
    print("Structure Modeling Complete!")
    print("=" * 80)
    print()
    print(f"Modeled structures saved in: {structures_dir}")
    print()
    print("Next steps:")
    print("  1. Inspect structures visually:")
    print(f"     pymol {structures_dir}/R175H/*.pdb")
    print()
    print("  2. Regenerate animations with real structures:")
    print("     python scripts/create_rescue_animation.py --use-modeled")
    print()
    print("  3. Calculate RMSD to quantify structural changes:")
    print("     pymol -c scripts/calculate_rmsd.pml")

    return 0


if __name__ == "__main__":
    sys.exit(main())
