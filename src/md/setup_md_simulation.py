#!/usr/bin/env python3
"""
Automate GROMACS MD simulation setup for p53 rescue mutants.

Usage:
    python src/md/setup_md_simulation.py --rescue "A189S,M133L,S95T" --target R175H

Requires:
    - GROMACS (gmx command)
    - Input PDB structure
"""

import argparse
import subprocess
import shutil
from pathlib import Path
from typing import List, Tuple


class MDSimulationSetup:
    """Automate MD simulation setup workflow."""

    def __init__(self, rescue_name: str, target: str, base_dir: Path):
        self.rescue_name = rescue_name.replace(",", "_")
        self.target = target
        self.base_dir = base_dir
        self.work_dir = base_dir / self.rescue_name
        self.work_dir.mkdir(parents=True, exist_ok=True)

        print(f"🧬 Setting up MD simulation for: {rescue_name}")
        print(f"   Target: {target}")
        print(f"   Work directory: {self.work_dir}")

    def find_input_structure(self) -> Path:
        """Locate input PDB structure for the rescue mutant."""

        # Try multiple possible locations
        search_paths = [
            Path(f"Data/processed/rescues/{self.target}/mutant_models/{self.rescue_name}.pdb"),
            Path(f"Data/processed/md_simulations/structures/{self.rescue_name}.pdb"),
            Path(f"Data/processed/cache/evoef2_repair/{self.rescue_name}_repaired.pdb"),
        ]

        for path in search_paths:
            if path.exists():
                print(f"  ✅ Found input structure: {path}")
                return path

        print(f"  ❌ Input structure not found. Searched:")
        for path in search_paths:
            print(f"     - {path}")
        raise FileNotFoundError(f"No input PDB found for {self.rescue_name}")

    def run_gmx_command(self, command: List[str], stdin: str = None) -> subprocess.CompletedProcess:
        """
        Run a GROMACS command.

        Args:
            command: List of command arguments
            stdin: Optional stdin input (for interactive commands)

        Returns:
            CompletedProcess object
        """
        try:
            print(f"  🔄 Running: {' '.join(command)}")
            result = subprocess.run(
                command,
                input=stdin,
                capture_output=True,
                text=True,
                check=True,
                cwd=self.work_dir
            )
            return result
        except subprocess.CalledProcessError as e:
            print(f"  ❌ Command failed: {e.cmd}")
            print(f"     STDOUT: {e.stdout}")
            print(f"     STDERR: {e.stderr}")
            raise

    def step1_pdb2gmx(self, input_pdb: Path):
        """
        Step 1: Generate topology using pdb2gmx.

        This step:
        - Adds hydrogen atoms
        - Generates topology file
        - Selects force field and water model
        """
        print("\n📋 Step 1: Generating topology (pdb2gmx)...")

        # Copy input PDB to work directory
        local_pdb = self.work_dir / "input.pdb"
        shutil.copy(input_pdb, local_pdb)

        command = [
            "gmx", "pdb2gmx",
            "-f", "input.pdb",
            "-o", "processed.gro",
            "-p", "topol.top",
            "-water", "spce",  # Use SPC/E instead of TIP3P for better compatibility
            "-ignh"  # Ignore existing hydrogens
        ]

        # Force field selection: 14 = AMBER99SB-ILDN
        # Water model selection: 1 = SPC/E
        stdin = "14\n1\n"

        try:
            self.run_gmx_command(command, stdin=stdin)
            print("  ✅ Topology generated: topol.top")
        except Exception as e:
            print(f"  ⚠️  Automated selection failed, trying with defaults...")
            # Try without stdin (use defaults)
            self.run_gmx_command(command[:-1])  # Remove -ignh

    def step2_editconf(self):
        """
        Step 2: Define simulation box.

        This step:
        - Creates cubic box around protein
        - Centers protein in box
        - Adds 1.0 nm padding
        """
        print("\n📦 Step 2: Defining simulation box (editconf)...")

        command = [
            "gmx", "editconf",
            "-f", "processed.gro",
            "-o", "box.gro",
            "-c",  # Center in box
            "-d", "1.0",  # Distance to box edge (nm)
            "-bt", "cubic"  # Cubic box
        ]

        self.run_gmx_command(command)
        print("  ✅ Simulation box defined: box.gro")

    def step3_solvate(self):
        """
        Step 3: Solvate with water.

        This step:
        - Fills box with TIP3P water molecules
        - Updates topology with water count
        """
        print("\n💧 Step 3: Solvating system...")

        command = [
            "gmx", "solvate",
            "-cp", "box.gro",
            "-cs", "spc216.gro",  # Water model
            "-o", "solvated.gro",
            "-p", "topol.top"
        ]

        self.run_gmx_command(command)
        print("  ✅ System solvated: solvated.gro")

    def step4_add_ions(self):
        """
        Step 4: Add ions for charge neutralization.

        This step:
        - Adds Na+ and Cl- ions
        - Neutralizes system charge
        - Adds physiological salt (0.15 M NaCl)
        """
        print("\n⚡ Step 4: Adding ions...")

        # First, generate .tpr file for genion
        # Use absolute path for MDP file
        project_root = Path.cwd()
        ions_mdp = project_root / "configs/md/mdp/ions.mdp"

        # Ensure MDP file exists
        if not ions_mdp.exists():
            ions_mdp.parent.mkdir(parents=True, exist_ok=True)
            ions_mdp.write_text("; Minimal MDP for ion generation\n")

        command_grompp = [
            "gmx", "grompp",
            "-f", str(ions_mdp),
            "-c", "solvated.gro",
            "-p", "topol.top",
            "-o", "ions.tpr",
            "-maxwarn", "1"
        ]

        self.run_gmx_command(command_grompp)

        # Now add ions
        command_genion = [
            "gmx", "genion",
            "-s", "ions.tpr",
            "-o", "system.gro",
            "-p", "topol.top",
            "-pname", "NA",
            "-nname", "CL",
            "-neutral",  # Neutralize charge
            "-conc", "0.15"  # 0.15 M NaCl
        ]

        # Select SOL (water) for ion replacement: group 13
        stdin = "13\n"
        self.run_gmx_command(command_genion, stdin=stdin)
        print("  ✅ Ions added: system.gro")

    def step5_energy_minimization(self):
        """
        Step 5: Energy minimization.

        This step:
        - Removes steric clashes
        - Optimizes geometry
        - Target: Fmax < 100 kJ/mol/nm
        """
        print("\n⚡ Step 5: Energy minimization...")

        # Use absolute path for MDP file
        project_root = Path.cwd()
        em_mdp = project_root / "configs/md/mdp/em.mdp"

        # Ensure MDP file exists
        if not em_mdp.exists():
            print("  ⚠️  em.mdp not found, creating default...")
            em_mdp.parent.mkdir(parents=True, exist_ok=True)
            em_mdp.write_text(self._default_em_mdp())

        # Generate .tpr
        command_grompp = [
            "gmx", "grompp",
            "-f", str(em_mdp),
            "-c", "system.gro",
            "-p", "topol.top",
            "-o", "em.tpr",
            "-maxwarn", "2"  # Allow warnings
        ]

        self.run_gmx_command(command_grompp)

        # Run energy minimization
        command_mdrun = [
            "gmx", "mdrun",
            "-v",  # Verbose
            "-deffnm", "em"
        ]

        self.run_gmx_command(command_mdrun)
        print("  ✅ Energy minimization complete: em.gro")

    def step6_equilibration_nvt(self):
        """
        Step 6: NVT equilibration (constant volume).

        This step:
        - Heats system to 300 K
        - 100 ps duration
        - Restraints on protein heavy atoms
        """
        print("\n🌡️  Step 6: NVT equilibration (100 ps)...")

        # Use absolute path for MDP file
        project_root = Path.cwd()
        nvt_mdp = project_root / "configs/md/mdp/nvt.mdp"

        # Ensure MDP file exists
        if not nvt_mdp.exists():
            print("  ⚠️  nvt.mdp not found, creating default...")
            nvt_mdp.parent.mkdir(parents=True, exist_ok=True)
            nvt_mdp.write_text(self._default_nvt_mdp())

        # Generate .tpr
        command_grompp = [
            "gmx", "grompp",
            "-f", str(nvt_mdp),
            "-c", "em.gro",
            "-r", "em.gro",  # Reference for restraints
            "-p", "topol.top",
            "-o", "nvt.tpr",
            "-maxwarn", "2"  # Allow warnings
        ]

        self.run_gmx_command(command_grompp)

        # Run NVT equilibration
        command_mdrun = [
            "gmx", "mdrun",
            "-v",
            "-deffnm", "nvt"
        ]

        self.run_gmx_command(command_mdrun)
        print("  ✅ NVT equilibration complete: nvt.gro")

    def step7_equilibration_npt(self):
        """
        Step 7: NPT equilibration (constant pressure).

        This step:
        - Equilibrates pressure to 1 bar
        - 100 ps duration
        - Maintains 300 K
        """
        print("\n🎚️  Step 7: NPT equilibration (100 ps)...")

        # Use absolute path for MDP file
        project_root = Path.cwd()
        npt_mdp = project_root / "configs/md/mdp/npt.mdp"

        # Ensure MDP file exists
        if not npt_mdp.exists():
            print("  ⚠️  npt.mdp not found, creating default...")
            npt_mdp.parent.mkdir(parents=True, exist_ok=True)
            npt_mdp.write_text(self._default_npt_mdp())

        # Generate .tpr
        command_grompp = [
            "gmx", "grompp",
            "-f", str(npt_mdp),
            "-c", "nvt.gro",
            "-r", "nvt.gro",
            "-t", "nvt.cpt",  # Continue from NVT
            "-p", "topol.top",
            "-o", "npt.tpr",
            "-maxwarn", "2"  # Allow warnings
        ]

        self.run_gmx_command(command_grompp)

        # Run NPT equilibration
        command_mdrun = [
            "gmx", "mdrun",
            "-v",
            "-deffnm", "npt"
        ]

        self.run_gmx_command(command_mdrun)
        print("  ✅ NPT equilibration complete: npt.gro")

    def step8_production_md(self, duration_ns: int = 50):
        """
        Step 8: Production MD simulation.

        Args:
            duration_ns: Simulation duration in nanoseconds (default: 50 ns)

        This step:
        - Runs full MD simulation
        - Saves trajectory every 10 ps
        - No restraints
        """
        print(f"\n🚀 Step 8: Production MD ({duration_ns} ns)...")
        print(f"   This will take ~{duration_ns / 4:.1f} hours on GPU (4 ns/hour)")
        print(f"   or ~{duration_ns * 2:.0f} hours on CPU (0.5 ns/hour)")

        # Use absolute path for MDP file
        project_root = Path.cwd()
        md_mdp = project_root / "configs/md/mdp/md.mdp"

        # Ensure MDP file exists
        if not md_mdp.exists():
            print("  ⚠️  md.mdp not found, creating default...")
            md_mdp.parent.mkdir(parents=True, exist_ok=True)
            md_mdp.write_text(self._default_md_mdp(duration_ns))

        # Generate .tpr
        command_grompp = [
            "gmx", "grompp",
            "-f", str(md_mdp),
            "-c", "npt.gro",
            "-t", "npt.cpt",
            "-p", "topol.top",
            "-o", "md.tpr",
            "-maxwarn", "2"  # Allow warnings
        ]

        self.run_gmx_command(command_grompp)

        # Run production MD
        print(f"\n  ⏳ Starting {duration_ns} ns MD simulation...")
        print(f"     Press Ctrl+C to abort (progress will be saved)")

        command_mdrun = [
            "gmx", "mdrun",
            "-v",
            "-deffnm", "md"
        ]

        try:
            self.run_gmx_command(command_mdrun)
            print(f"  ✅ Production MD complete: md.xtc ({duration_ns} ns)")
        except KeyboardInterrupt:
            print("\n  ⚠️  Simulation interrupted by user")
            print("     Restart with: gmx mdrun -v -deffnm md -cpi md.cpt")

    @staticmethod
    def _default_em_mdp() -> str:
        """Default energy minimization MDP parameters."""
        return """
; Energy minimization
integrator  = steep         ; Steepest descent
emtol       = 100.0         ; Stop when Fmax < 100 kJ/mol/nm
emstep      = 0.01          ; Initial step size
nsteps      = 5000          ; Maximum steps

; Output control
nstxout     = 100
nstvout     = 100
nstfout     = 100
nstlog      = 100
nstenergy   = 100

; Neighbor searching
cutoff-scheme   = Verlet
ns_type         = grid
nstlist         = 10
rlist           = 1.0

; Electrostatics
coulombtype     = PME
rcoulomb        = 1.0

; VdW
vdwtype         = cutoff
rvdw            = 1.0

; PME
pme_order       = 4
fourierspacing  = 0.16
"""

    @staticmethod
    def _default_nvt_mdp() -> str:
        """Default NVT equilibration MDP parameters."""
        return """
; NVT equilibration
define      = -DPOSRES      ; Position restraints
integrator  = md            ; Leap-frog
dt          = 0.002         ; 2 fs timestep
nsteps      = 50000         ; 100 ps

; Output
nstxout     = 500
nstvout     = 500
nstenergy   = 500
nstlog      = 500

; Neighbor searching
cutoff-scheme   = Verlet
ns_type         = grid
nstlist         = 10
rlist           = 1.0

; Electrostatics
coulombtype     = PME
rcoulomb        = 1.0

; VdW
vdwtype         = cutoff
rvdw            = 1.0

; Temperature coupling
tcoupl          = V-rescale
tc-grps         = Protein Non-Protein
tau_t           = 0.1 0.1
ref_t           = 300 300

; Pressure coupling
pcoupl          = no

; Velocity generation
gen_vel         = yes
gen_temp        = 300
gen_seed        = -1

; Constraints
constraints     = h-bonds
constraint_algorithm = lincs
"""

    @staticmethod
    def _default_npt_mdp() -> str:
        """Default NPT equilibration MDP parameters."""
        return """
; NPT equilibration
define      = -DPOSRES
integrator  = md
dt          = 0.002
nsteps      = 50000         ; 100 ps

; Output
nstxout     = 500
nstvout     = 500
nstenergy   = 500
nstlog      = 500

; Neighbor searching
cutoff-scheme   = Verlet
ns_type         = grid
nstlist         = 10
rlist           = 1.0

; Electrostatics
coulombtype     = PME
rcoulomb        = 1.0

; VdW
vdwtype         = cutoff
rvdw            = 1.0

; Temperature coupling
tcoupl          = V-rescale
tc-grps         = Protein Non-Protein
tau_t           = 0.1 0.1
ref_t           = 300 300

; Pressure coupling
pcoupl              = Parrinello-Rahman
pcoupltype          = isotropic
tau_p               = 2.0
ref_p               = 1.0
compressibility     = 4.5e-5

; Velocity generation
gen_vel         = no

; Constraints
constraints     = h-bonds
constraint_algorithm = lincs
"""

    @staticmethod
    def _default_md_mdp(duration_ns: int) -> str:
        """Default production MD MDP parameters."""
        nsteps = duration_ns * 500000  # 2 fs timestep, 500000 steps = 1 ns
        return f"""
; Production MD ({duration_ns} ns)
integrator  = md
dt          = 0.002
nsteps      = {nsteps}

; Output (save every 10 ps)
nstxout     = 0             ; Don't save coordinates (use xtc)
nstvout     = 0
nstxtcout   = 5000          ; Save to xtc every 10 ps
nstenergy   = 5000
nstlog      = 5000

; Neighbor searching
cutoff-scheme   = Verlet
ns_type         = grid
nstlist         = 10
rlist           = 1.0

; Electrostatics
coulombtype     = PME
rcoulomb        = 1.0

; VdW
vdwtype         = cutoff
rvdw            = 1.0

; Temperature coupling
tcoupl          = V-rescale
tc-grps         = Protein Non-Protein
tau_t           = 0.1 0.1
ref_t           = 300 300

; Pressure coupling
pcoupl              = Parrinello-Rahman
pcoupltype          = isotropic
tau_p               = 2.0
ref_p               = 1.0
compressibility     = 4.5e-5

; Velocity generation
gen_vel         = no

; Constraints
constraints     = h-bonds
constraint_algorithm = lincs

; PBC
pbc             = xyz
"""

    def run_full_pipeline(self, duration_ns: int = 50):
        """Run complete MD simulation pipeline."""
        try:
            # Find input structure
            input_pdb = self.find_input_structure()

            # Run all steps
            self.step1_pdb2gmx(input_pdb)
            self.step2_editconf()
            self.step3_solvate()
            self.step4_add_ions()
            self.step5_energy_minimization()
            self.step6_equilibration_nvt()
            self.step7_equilibration_npt()
            self.step8_production_md(duration_ns)

            print("\n" + "=" * 80)
            print("✅ MD Simulation Setup Complete!")
            print("=" * 80)
            print(f"Output directory: {self.work_dir}")
            print(f"Trajectory: {self.work_dir / 'md.xtc'}")
            print(f"\nNext steps:")
            print(f"1. Analyze trajectory: python src/md/analyze_trajectory.py")
            print(f"2. Extract conformers for docking")

        except Exception as e:
            print(f"\n❌ Pipeline failed: {e}")
            raise


def main():
    parser = argparse.ArgumentParser(description="Setup GROMACS MD simulation for p53 rescue mutant")
    parser.add_argument("--rescue", required=True, help="Rescue mutation (e.g., 'A189S,M133L,S95T')")
    parser.add_argument("--target", required=True, help="Target mutation (e.g., 'R175H')")
    parser.add_argument("--duration", type=int, default=50, help="Simulation duration in ns (default: 50)")
    parser.add_argument("--base-dir", type=Path, default=Path("Data/processed/md_simulations"),
                        help="Base directory for MD simulations")

    args = parser.parse_args()

    # Create setup instance
    setup = MDSimulationSetup(args.rescue, args.target, args.base_dir)

    # Run pipeline
    setup.run_full_pipeline(duration_ns=args.duration)


if __name__ == "__main__":
    main()
