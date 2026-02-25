#!/usr/bin/env python
"""
p53-proteoMgCAD MD Validation Script
Generated: 2026-02-25T07:03:09.440819
Target: R248W_R249S
Mutations: L111V, C124W, A138N, C141L, L145F, V157I, V203A, N235K, Y236F, N239Q, C242S, R248W, R249S, T253V, I255F, G262E, N268R, V272A, V274I, C275S, A276S
"""

from openmm import *
from openmm.app import *
from openmm.unit import *
import numpy as np
import sys

# ============================================================================
# CONFIGURATION
# ============================================================================

SEQUENCE = """MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRVGFLHSGTAKSVTWTYSPALNKMFCQLNKTLPVQFWVDSTPPPGTRIRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRAEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYKFMCQSSSMGGMNWSPILVIFTLEDSSENLLGRRSFEARISSCPGRDRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPGGSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD"""
NAME = "R248W_R249S"
MUTATIONS = ['L111V', 'C124W', 'A138N', 'C141L', 'L145F', 'V157I', 'V203A', 'N235K', 'Y236F', 'N239Q', 'C242S', 'R248W', 'R249S', 'T253V', 'I255F', 'G262E', 'N268R', 'V272A', 'V274I', 'C275S', 'A276S']

# Simulation parameters
TEMPERATURE = 310.0 * kelvin
PRESSURE = 1.0 * bar
TIMESTEP = 0.002 * picoseconds
EQUILIBRATION_STEPS = 50000
PRODUCTION_STEPS = 5000000
OUTPUT_INTERVAL = 5000

# ============================================================================
# SETUP
# ============================================================================

print(f"Setting up MD simulation for {NAME}")
print(f"Sequence length: {len(SEQUENCE)} residues")
print(f"Mutations: {MUTATIONS}")

# Load structure (from ESMFold prediction or PDB)
# In practice, you would load from file:
# pdb = PDBFile('input.pdb')

# For this template, we assume structure is provided
try:
    pdb = PDBFile(f'{NAME}.pdb')
except:
    print("ERROR: Please provide {NAME}.pdb file")
    print("You can generate this using ESMFold or AlphaFold")
    sys.exit(1)

# Force field
forcefield = ForceField('amber14-all.xml', 'tip3p.xml')

# Create modeller for solvation
modeller = Modeller(pdb.topology, pdb.positions)

# Add hydrogens
modeller.addHydrogens(forcefield)

# Add solvent
print("Adding solvent box...")
modeller.addSolvent(
    forcefield,
    model='tip3p',
    padding=1.0*nanometers,
    ionicStrength=0.15*molar
)

# Create system
print("Creating system...")
system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=1.0*nanometers,
    constraints=HBonds
)

# Add barostat for NPT
system.addForce(MonteCarloBarostat(PRESSURE, TEMPERATURE))

# Integrator
integrator = LangevinMiddleIntegrator(TEMPERATURE, 1/picosecond, TIMESTEP)

# Platform selection
platform_name = "auto"
if platform_name == "auto":
    try:
        platform = Platform.getPlatformByName('CUDA')
    except:
        try:
            platform = Platform.getPlatformByName('OpenCL')
        except:
            platform = Platform.getPlatformByName('CPU')
else:
    platform = Platform.getPlatformByName(platform_name)

print(f"Using platform: {platform.getName()}")

# Create simulation
simulation = Simulation(modeller.topology, system, integrator, platform)
simulation.context.setPositions(modeller.positions)

# ============================================================================
# MINIMIZATION
# ============================================================================

print("Minimizing energy...")
simulation.minimizeEnergy(maxIterations=1000)

# ============================================================================
# EQUILIBRATION
# ============================================================================

print(f"Equilibrating for {EQUILIBRATION_STEPS} steps...")
simulation.context.setVelocitiesToTemperature(TEMPERATURE)

# NVT equilibration (no barostat)
simulation.step(EQUILIBRATION_STEPS // 2)

# NPT equilibration
simulation.step(EQUILIBRATION_STEPS // 2)

# ============================================================================
# PRODUCTION
# ============================================================================

print(f"Running production MD for {PRODUCTION_STEPS} steps...")
print(f"Output interval: {OUTPUT_INTERVAL} steps")

# Reporters
simulation.reporters.append(
    DCDReporter(f'{NAME}_trajectory.dcd', OUTPUT_INTERVAL)
)
simulation.reporters.append(
    StateDataReporter(
        f'{NAME}_energy.csv',
        OUTPUT_INTERVAL,
        step=True,
        time=True,
        potentialEnergy=True,
        kineticEnergy=True,
        temperature=True,
        volume=True,
        density=True
    )
)
simulation.reporters.append(
    StateDataReporter(sys.stdout, OUTPUT_INTERVAL * 10, step=True, time=True,
                     remainingTime=True, speed=True, totalSteps=PRODUCTION_STEPS)
)

# Run production
simulation.step(PRODUCTION_STEPS)

print("\nSimulation complete!")
print(f"Trajectory: {NAME}_trajectory.dcd")
print(f"Energy: {NAME}_energy.csv")

# ============================================================================
# ANALYSIS (Basic)
# ============================================================================

print("\nRunning basic analysis...")

import mdtraj as md

# Load trajectory
traj = md.load(f'{NAME}_trajectory.dcd', top=f'{NAME}.pdb')

# RMSD
rmsd = md.rmsd(traj, traj[0], atom_indices=traj.topology.select('name CA'))
print(f"RMSD: {rmsd.mean()*10:.2f} ± {rmsd.std()*10:.2f} Å")

# RMSF
rmsf = md.rmsf(traj, traj[0], atom_indices=traj.topology.select('name CA'))
print(f"Mean RMSF: {rmsf.mean()*10:.2f} Å")

# Radius of gyration
rg = md.compute_rg(traj)
print(f"Radius of gyration: {rg.mean()*10:.2f} ± {rg.std()*10:.2f} Å")

# Save analysis
analysis = {
    "name": NAME,
    "mutations": MUTATIONS,
    "rmsd_mean": float(rmsd.mean()),
    "rmsd_std": float(rmsd.std()),
    "rmsf_mean": float(rmsf.mean()),
    "rg_mean": float(rg.mean()),
    "rg_std": float(rg.std()),
    "n_frames": len(traj),
    "simulation_time_ns": len(traj) * TIMESTEP.value_in_unit(picoseconds) * OUTPUT_INTERVAL / 1000
}

import json
with open(f'{NAME}_analysis.json', 'w') as f:
    json.dump(analysis, f, indent=2)

print(f"\nAnalysis saved to {NAME}_analysis.json")
