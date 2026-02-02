#!/usr/bin/env python3
"""
p53 Explicit Solvent MD Simulation Script

Run locally with:
    pip install openmm pdbfixer mdtraj numpy matplotlib
    python scripts/run_md_explicit.py

For GPU acceleration, install CUDA toolkit and openmm-cuda.
"""
from __future__ import annotations

import os
import sys
import time
import json
import argparse
from pathlib import Path

import numpy as np

# Check for required packages
try:
    from openmm import *
    from openmm.app import *
    from openmm.unit import *
except ImportError:
    print("ERROR: OpenMM not installed. Run: pip install openmm")
    sys.exit(1)

try:
    from pdbfixer import PDBFixer
except ImportError:
    print("ERROR: PDBFixer not installed. Run: pip install pdbfixer")
    sys.exit(1)

try:
    import mdtraj as md
except ImportError:
    print("ERROR: MDTraj not installed. Run: pip install mdtraj")
    sys.exit(1)


# ============================================================================
# Configuration
# ============================================================================

# p53 sequence
P53_FULL = (
    "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP"
    "DEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAK"
    "SVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHE"
    "RCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNS"
    "SCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELP"
    "PGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPG"
    "GSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD"
)
P53_CORE = P53_FULL[93:312]  # Core domain (residues 94-312)
CORE_START = 94

# Default simulation parameters
DEFAULT_CONFIG = {
    "timestep_fs": 2.0,
    "equilibration_ps": 100,
    "production_ns": 1.0,
    "save_interval_ps": 2.0,
    "temperature_k": 300,
    "pressure_bar": 1.0,
    "padding_nm": 1.0,
    "ionic_strength_M": 0.15,
    "platform": "CPU",  # Options: CPU, OpenCL, CUDA
}

# Variants to simulate
DEFAULT_VARIANTS = {
    "WT": [],
    "R175H": ["R175H"],
    "R175H_N239Y": ["R175H", "N239Y"],
}


# ============================================================================
# Helper Functions
# ============================================================================

def apply_mutations(sequence: str, mutations: list[str]) -> str:
    """Apply mutations to the p53 core domain sequence."""
    for mut in mutations:
        wt_aa = mut[0]
        pos = int(mut[1:-1])
        mut_aa = mut[-1]
        core_pos = pos - CORE_START

        if core_pos < 0 or core_pos >= len(sequence):
            raise ValueError(f"Position {pos} out of range for core domain")
        if sequence[core_pos] != wt_aa:
            raise ValueError(f"Expected {wt_aa} at position {pos}, found {sequence[core_pos]}")

        sequence = sequence[:core_pos] + mut_aa + sequence[core_pos + 1:]
        print(f"  Applied {mut} at core position {core_pos}")

    return sequence


def get_structure_from_esmfold(sequence: str, output_path: Path, max_retries: int = 3) -> Path:
    """Predict structure using ESMFold API."""
    import requests

    url = "https://api.esmatlas.com/foldSequence/v1/pdb/"

    for attempt in range(max_retries):
        try:
            print(f"  Calling ESMFold API (attempt {attempt + 1}/{max_retries})...")
            response = requests.post(
                url,
                data=sequence,
                headers={"Content-Type": "text/plain"},
                timeout=300,
            )
            if response.status_code == 200:
                output_path.write_text(response.text)
                print(f"  Success! Saved to {output_path}")
                return output_path
            elif response.status_code == 503:
                print(f"  Server busy, waiting 30s...")
                time.sleep(30)
            else:
                print(f"  Error: HTTP {response.status_code}")
                time.sleep(10)
        except requests.Timeout:
            print(f"  Timeout, retrying...")
            time.sleep(10)
        except Exception as e:
            print(f"  Error: {e}")
            time.sleep(10)

    raise RuntimeError("ESMFold API failed after all retries")


def fix_structure(pdb_path: Path) -> PDBFixer:
    """Fix structure with PDBFixer (add missing atoms, hydrogens)."""
    fixer = PDBFixer(filename=str(pdb_path))
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    return fixer


def select_platform(preferred: str) -> Platform:
    """Select OpenMM platform."""
    available = [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())]
    print(f"  Available platforms: {available}")

    if preferred in available:
        platform = Platform.getPlatformByName(preferred)
        print(f"  Using: {preferred}")
        return platform

    # Fallback order
    for name in ["CUDA", "OpenCL", "CPU", "Reference"]:
        if name in available:
            platform = Platform.getPlatformByName(name)
            print(f"  Using fallback: {name}")
            return platform

    raise RuntimeError("No suitable platform found")


# ============================================================================
# Main Simulation Function
# ============================================================================

def run_simulation(
    name: str,
    mutations: list[str],
    config: dict,
    output_dir: Path,
    use_existing_structure: Path | None = None,
) -> dict:
    """
    Run explicit solvent MD simulation.

    Args:
        name: Variant name (e.g., 'WT', 'R175H')
        mutations: List of mutations to apply
        config: Simulation configuration
        output_dir: Output directory
        use_existing_structure: Path to existing PDB (skip ESMFold if provided)

    Returns:
        Dictionary with simulation results
    """
    print(f"\n{'=' * 60}")
    print(f"SIMULATING: {name}")
    print(f"{'=' * 60}")

    struct_dir = output_dir / "structures"
    traj_dir = output_dir / "trajectories"
    struct_dir.mkdir(parents=True, exist_ok=True)
    traj_dir.mkdir(parents=True, exist_ok=True)

    # Convert config to simulation units
    timestep = config["timestep_fs"] * femtoseconds
    eq_steps = int(config["equilibration_ps"] * 1000 / config["timestep_fs"])
    prod_steps = int(config["production_ns"] * 1e6 / config["timestep_fs"])
    save_interval = int(config["save_interval_ps"] * 1000 / config["timestep_fs"])
    temperature = config["temperature_k"] * kelvin
    pressure = config["pressure_bar"] * bar

    # 1. Prepare sequence
    print(f"\n[{name}] Preparing sequence...")
    sequence = P53_CORE
    if mutations:
        sequence = apply_mutations(sequence, mutations)
    else:
        print("  No mutations (wild-type)")

    # 2. Get/predict structure
    raw_pdb = struct_dir / f"{name}_raw.pdb"
    if use_existing_structure and use_existing_structure.exists():
        print(f"\n[{name}] Using existing structure: {use_existing_structure}")
        import shutil
        shutil.copy(use_existing_structure, raw_pdb)
    elif raw_pdb.exists():
        print(f"\n[{name}] Using cached structure: {raw_pdb}")
    else:
        print(f"\n[{name}] Predicting structure with ESMFold...")
        get_structure_from_esmfold(sequence, raw_pdb)

    # 3. Fix structure
    print(f"\n[{name}] Fixing structure with PDBFixer...")
    fixer = fix_structure(raw_pdb)
    print(f"  Fixed structure: {fixer.topology.getNumAtoms()} atoms")

    # 4. Create explicit solvent system
    print(f"\n[{name}] Creating explicit solvent system...")
    forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")

    modeller = Modeller(fixer.topology, fixer.positions)
    print(f"  Adding water box ({config['padding_nm']} nm padding)...")
    modeller.addSolvent(
        forcefield,
        model="tip3p",
        padding=config["padding_nm"] * nanometer,
        ionicStrength=config["ionic_strength_M"] * molar,
    )
    n_atoms = modeller.topology.getNumAtoms()
    print(f"  Solvated system: {n_atoms} atoms")

    # Save solvated structure
    solvated_pdb = struct_dir / f"{name}_solvated.pdb"
    with open(solvated_pdb, "w") as f:
        PDBFile.writeFile(modeller.topology, modeller.positions, f)

    # Create system
    print(f"  Creating OpenMM system...")
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1.0 * nanometer,
        constraints=HBonds,
    )

    # Select platform
    platform = select_platform(config["platform"])

    # 5. Energy minimization
    print(f"\n[{name}] Energy minimization...")
    integrator = LangevinMiddleIntegrator(temperature, 1 / picosecond, timestep)
    simulation = Simulation(modeller.topology, system, integrator, platform)
    simulation.context.setPositions(modeller.positions)

    e_initial = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    simulation.minimizeEnergy(maxIterations=1000)
    e_final = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    print(f"  Energy: {e_initial} -> {e_final}")

    positions = simulation.context.getState(getPositions=True).getPositions()

    # 6. NPT Equilibration
    print(f"\n[{name}] NPT Equilibration ({config['equilibration_ps']} ps)...")

    # Add barostat
    system.addForce(MonteCarloBarostat(pressure, temperature))

    integrator = LangevinMiddleIntegrator(temperature, 1 / picosecond, timestep)
    simulation = Simulation(modeller.topology, system, integrator, platform)
    simulation.context.setPositions(positions)
    simulation.context.setVelocitiesToTemperature(temperature)

    simulation.reporters.append(
        StateDataReporter(
            sys.stdout, eq_steps // 5,
            step=True, temperature=True, progress=True,
            remainingTime=True, speed=True, totalSteps=eq_steps,
        )
    )

    simulation.step(eq_steps)

    # Save equilibrated structure
    eq_positions = simulation.context.getState(getPositions=True).getPositions()
    eq_pdb = struct_dir / f"{name}_equilibrated.pdb"
    with open(eq_pdb, "w") as f:
        PDBFile.writeFile(modeller.topology, eq_positions, f)
    print(f"  Saved: {eq_pdb}")

    # 7. Production MD
    print(f"\n[{name}] Production MD ({config['production_ns']} ns)...")

    simulation.reporters.clear()

    traj_path = traj_dir / f"{name}_traj.dcd"
    log_path = traj_dir / f"{name}_log.csv"

    if traj_path.exists():
        traj_path.unlink()

    simulation.reporters.append(DCDReporter(str(traj_path), save_interval))
    simulation.reporters.append(
        StateDataReporter(
            str(log_path), save_interval,
            step=True, time=True, potentialEnergy=True, temperature=True,
        )
    )
    simulation.reporters.append(
        StateDataReporter(
            sys.stdout, prod_steps // 10,
            step=True, time=True, potentialEnergy=True, temperature=True,
            progress=True, remainingTime=True, speed=True, totalSteps=prod_steps,
        )
    )

    simulation.step(prod_steps)

    # 8. Analysis
    print(f"\n[{name}] Analyzing trajectory...")
    traj = md.load(str(traj_path), top=str(eq_pdb))

    # Select protein atoms only
    protein_idx = traj.topology.select("protein")
    protein_traj = traj.atom_slice(protein_idx)
    print(f"  Protein atoms: {protein_traj.n_atoms}")

    # RMSD
    rmsd = md.rmsd(protein_traj, protein_traj, 0) * 10  # nm to Angstrom

    # RMSF (CA atoms)
    protein_aligned = protein_traj.superpose(protein_traj, 0)
    ca_idx = protein_traj.topology.select("name CA")
    rmsf = md.rmsf(protein_aligned, protein_aligned, frame=0, atom_indices=ca_idx) * 10
    residue_nums = [protein_traj.topology.atom(i).residue.resSeq for i in ca_idx]

    # Save analysis
    np.save(traj_dir / f"{name}_rmsd.npy", rmsd)
    np.save(traj_dir / f"{name}_rmsf.npy", rmsf)

    results = {
        "name": name,
        "mutations": mutations,
        "n_frames": traj.n_frames,
        "n_atoms_protein": protein_traj.n_atoms,
        "n_atoms_total": n_atoms,
        "mean_rmsd": float(np.mean(rmsd)),
        "final_rmsd": float(rmsd[-1]),
        "max_rmsd": float(np.max(rmsd)),
        "mean_rmsf": float(np.mean(rmsf)),
    }

    print(f"\n[{name}] Complete!")
    print(f"  Frames: {traj.n_frames}")
    print(f"  Mean RMSD: {results['mean_rmsd']:.2f} A")
    print(f"  Final RMSD: {results['final_rmsd']:.2f} A")

    return results


# ============================================================================
# Main Entry Point
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Run explicit solvent MD simulation for p53 variants"
    )
    parser.add_argument(
        "--variants", nargs="+", default=["WT", "R175H", "R175H_N239Y"],
        help="Variants to simulate (default: WT R175H R175H_N239Y)"
    )
    parser.add_argument(
        "--production-ns", type=float, default=1.0,
        help="Production simulation time in ns (default: 1.0)"
    )
    parser.add_argument(
        "--equilibration-ps", type=float, default=100,
        help="Equilibration time in ps (default: 100)"
    )
    parser.add_argument(
        "--platform", choices=["CPU", "OpenCL", "CUDA"], default="CPU",
        help="OpenMM platform (default: CPU)"
    )
    parser.add_argument(
        "--output-dir", type=Path, default=Path("Data/processed/md_simulations"),
        help="Output directory"
    )
    parser.add_argument(
        "--use-alphafold", type=Path, default=None,
        help="Use existing AlphaFold structure instead of ESMFold"
    )

    args = parser.parse_args()

    # Build config
    config = DEFAULT_CONFIG.copy()
    config["production_ns"] = args.production_ns
    config["equilibration_ps"] = args.equilibration_ps
    config["platform"] = args.platform

    # Print header
    print("=" * 60)
    print("p53 EXPLICIT SOLVENT MD SIMULATION")
    print("=" * 60)
    print(f"Variants: {args.variants}")
    print(f"Production: {config['production_ns']} ns")
    print(f"Equilibration: {config['equilibration_ps']} ps")
    print(f"Platform: {config['platform']}")
    print(f"Output: {args.output_dir}")
    print("=" * 60)

    # Run simulations
    all_results = {}
    for name in args.variants:
        if name not in DEFAULT_VARIANTS:
            # Parse custom variant (e.g., "R248Q" or "R175H_A159V")
            mutations = name.split("_") if "_" in name else [name] if name != "WT" else []
        else:
            mutations = DEFAULT_VARIANTS[name]

        results = run_simulation(
            name=name,
            mutations=mutations,
            config=config,
            output_dir=args.output_dir,
            use_existing_structure=args.use_alphafold,
        )
        all_results[name] = results

    # Summary
    print("\n" + "=" * 60)
    print("SIMULATION SUMMARY")
    print("=" * 60)

    for name, r in all_results.items():
        status = "STABLE" if r["final_rmsd"] < 2.5 else "UNSTABLE"
        print(f"\n{name}:")
        print(f"  Final RMSD: {r['final_rmsd']:.2f} A")
        print(f"  Mean RMSD:  {r['mean_rmsd']:.2f} A")
        print(f"  Status:     {status}")

    # Save results
    results_path = args.output_dir / "md_results.json"
    with open(results_path, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\nResults saved to: {results_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
