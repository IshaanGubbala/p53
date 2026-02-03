#!/usr/bin/env python3
"""
Run ensemble docking for rescue-drug pairs.

Usage:
    python src/docking/run_ensemble_docking.py \
        --rescue "A189S_M133L_S95T" \
        --drug COTI-2

Requires:
    - AutoDock Vina
    - Prepared drug structures (PDBQT)
    - Representative conformers from MD analysis
"""

import argparse
import json
import subprocess
from pathlib import Path
from typing import List, Dict, Tuple
import numpy as np


class EnsembleDocking:
    """Automated ensemble docking with AutoDock Vina."""

    def __init__(self, rescue_name: str, drug_name: str):
        self.rescue_name = rescue_name
        self.drug_name = drug_name

        # Directories
        self.md_dir = Path("Data/processed/md_simulations") / rescue_name
        self.conformers_dir = self.md_dir / "analysis" / "conformers"
        self.drug_dir = Path("Data/raw/drugs")
        self.output_dir = Path("Data/processed/docking/results") / rescue_name / drug_name
        self.output_dir.mkdir(parents=True, exist_ok=True)

        print(f"🎯 Setting up ensemble docking")
        print(f"   Rescue: {rescue_name}")
        print(f"   Drug: {drug_name}")
        print(f"   Output: {self.output_dir}")

    def find_drug_structure(self) -> Path:
        """Locate drug structure file (SDF or MOL2 or PDBQT)."""

        # Check multiple locations and formats
        search_paths = [
            self.drug_dir / "prepared" / f"{self.drug_name}.pdbqt",
            self.drug_dir / "structures" / f"{self.drug_name}.pdb",
            self.drug_dir / "structures" / f"{self.drug_name}.sdf",
            self.drug_dir / "structures" / f"{self.drug_name}.mol2",
        ]

        for path in search_paths:
            if path.exists():
                print(f"   ✅ Found drug structure: {path}")
                return path

        raise FileNotFoundError(
            f"Drug structure not found: {self.drug_name}\n"
            f"Searched: {search_paths}"
        )

    def find_conformers(self) -> List[Path]:
        """Find all conformer PDB files."""

        if not self.conformers_dir.exists():
            raise FileNotFoundError(
                f"Conformers directory not found: {self.conformers_dir}\n"
                f"Run trajectory analysis first: python src/md/analyze_trajectory.py"
            )

        conformers = sorted(self.conformers_dir.glob("conformer_*.pdb"))

        if not conformers:
            raise FileNotFoundError(f"No conformers found in {self.conformers_dir}")

        print(f"   ✅ Found {len(conformers)} conformers")
        return conformers

    def prepare_receptor(self, pdb_path: Path) -> Path:
        """
        Prepare receptor PDB for docking (convert to PDBQT).

        Args:
            pdb_path: Path to conformer PDB file

        Returns:
            Path to PDBQT file
        """

        pdbqt_path = self.output_dir / f"{pdb_path.stem}.pdbqt"

        if pdbqt_path.exists():
            # Already prepared
            return pdbqt_path

        try:
            # Use prepare_receptor4.py from AutoDockTools
            cmd = [
                "prepare_receptor4.py",
                "-r", str(pdb_path),
                "-o", str(pdbqt_path),
                "-A", "hydrogens",
                "-U", "nphs"  # Merge non-polar hydrogens
            ]

            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )

            return pdbqt_path

        except subprocess.CalledProcessError:
            # Fallback: Try using obabel if prepare_receptor4.py not available
            print(f"   ⚠️  prepare_receptor4.py not available, using obabel...")
            try:
                cmd = [
                    "obabel",
                    str(pdb_path),
                    "-O", str(pdbqt_path),
                    "-xr"  # Rigid receptor
                ]
                subprocess.run(cmd, capture_output=True, text=True, check=True)
                return pdbqt_path
            except:
                raise RuntimeError(
                    f"Failed to prepare receptor: {pdb_path}\n"
                    f"Install AutoDockTools or check obabel installation"
                )

    def prepare_ligand(self, drug_path: Path) -> Path:
        """
        Prepare ligand for docking (convert to PDBQT if needed).

        Args:
            drug_path: Path to drug structure file

        Returns:
            Path to PDBQT file
        """

        # If already PDBQT, return as-is
        if drug_path.suffix == ".pdbqt":
            return drug_path

        # Otherwise, convert
        pdbqt_path = self.output_dir / f"{self.drug_name}.pdbqt"

        if pdbqt_path.exists():
            return pdbqt_path

        try:
            # Try prepare_ligand4.py first
            cmd = [
                "prepare_ligand4.py",
                "-l", str(drug_path),
                "-o", str(pdbqt_path),
                "-A", "hydrogens"
            ]

            subprocess.run(cmd, capture_output=True, text=True, check=True)
            return pdbqt_path

        except:
            # Fallback to obabel
            print(f"   ⚠️  prepare_ligand4.py not available, using obabel...")
            try:
                cmd = [
                    "obabel",
                    str(drug_path),
                    "-O", str(pdbqt_path),
                    "-h"  # Add hydrogens
                ]
                subprocess.run(cmd, capture_output=True, text=True, check=True)
                return pdbqt_path
            except Exception as e:
                raise RuntimeError(f"Failed to prepare ligand: {e}")

    def define_binding_site(self, conformer_pdb: Path) -> Dict[str, float]:
        """
        Define binding site coordinates.

        For now, uses a simple center-of-mass approach.
        In production, should use known binding sites or blind docking.

        Args:
            conformer_pdb: Path to conformer PDB

        Returns:
            Dict with center_x, center_y, center_z, size_x, size_y, size_z
        """

        # Parse PDB to get protein center of mass
        coords = []
        with open(conformer_pdb) as f:
            for line in f:
                if line.startswith("ATOM"):
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append([x, y, z])

        coords = np.array(coords)
        center = coords.mean(axis=0)

        # Define box size (20 Å cube around center)
        box_size = 20.0

        return {
            "center_x": float(center[0]),
            "center_y": float(center[1]),
            "center_z": float(center[2]),
            "size_x": box_size,
            "size_y": box_size,
            "size_z": box_size
        }

    def run_vina_docking(
        self,
        receptor_pdbqt: Path,
        ligand_pdbqt: Path,
        binding_site: Dict[str, float],
        conformer_id: int
    ) -> Tuple[float, Path]:
        """
        Run AutoDock Vina docking.

        Args:
            receptor_pdbqt: Receptor PDBQT file
            ligand_pdbqt: Ligand PDBQT file
            binding_site: Binding site coordinates
            conformer_id: Conformer number

        Returns:
            (best_affinity, output_pdbqt)
        """

        output_pdbqt = self.output_dir / f"conformer_{conformer_id}_docked.pdbqt"
        log_file = self.output_dir / f"conformer_{conformer_id}_log.txt"

        # Create config file
        config_file = self.output_dir / f"conformer_{conformer_id}_config.txt"
        config_content = f"""receptor = {receptor_pdbqt}
ligand = {ligand_pdbqt}
center_x = {binding_site['center_x']:.3f}
center_y = {binding_site['center_y']:.3f}
center_z = {binding_site['center_z']:.3f}
size_x = {binding_site['size_x']:.1f}
size_y = {binding_site['size_y']:.1f}
size_z = {binding_site['size_z']:.1f}
exhaustiveness = 32
num_modes = 10
energy_range = 3
"""
        config_file.write_text(config_content)

        # Run Vina
        cmd = [
            "vina",
            "--config", str(config_file),
            "--out", str(output_pdbqt),
            "--log", str(log_file)
        ]

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True,
                timeout=300  # 5 min timeout per docking
            )

            # Parse best affinity from log
            with open(log_file) as f:
                log_content = f.read()

            # Find best affinity (first line in results)
            for line in log_content.split('\n'):
                if line.strip().startswith('1 '):
                    parts = line.split()
                    best_affinity = float(parts[1])
                    return best_affinity, output_pdbqt

            raise RuntimeError("Could not parse affinity from Vina log")

        except subprocess.TimeoutExpired:
            print(f"   ⚠️  Docking timed out for conformer {conformer_id}")
            return None, None
        except Exception as e:
            print(f"   ❌ Docking failed for conformer {conformer_id}: {e}")
            return None, None

    def run_ensemble_docking(self) -> Dict:
        """
        Run docking on all conformers and calculate consensus.

        Returns:
            Dict with results
        """

        print(f"\n🚀 Starting ensemble docking...")

        # Find drug and conformers
        drug_path = self.find_drug_structure()
        conformers = self.find_conformers()

        # Prepare ligand once
        print(f"\n📋 Preparing ligand...")
        ligand_pdbqt = self.prepare_ligand(drug_path)
        print(f"   ✅ Ligand prepared: {ligand_pdbqt}")

        # Dock to each conformer
        affinities = []
        successful_dockings = 0

        for i, conformer_pdb in enumerate(conformers, start=1):
            print(f"\n🎯 Docking to conformer {i}/{len(conformers)}...")

            # Prepare receptor
            receptor_pdbqt = self.prepare_receptor(conformer_pdb)

            # Define binding site
            binding_site = self.define_binding_site(conformer_pdb)

            # Run docking
            affinity, output_pdbqt = self.run_vina_docking(
                receptor_pdbqt,
                ligand_pdbqt,
                binding_site,
                i
            )

            if affinity is not None:
                affinities.append(affinity)
                successful_dockings += 1
                print(f"   ✅ Affinity: {affinity:.2f} kcal/mol")
            else:
                print(f"   ❌ Docking failed")

        # Calculate consensus
        if not affinities:
            raise RuntimeError("All dockings failed!")

        affinities = np.array(affinities)
        consensus_affinity = np.median(affinities)
        affinity_std = np.std(affinities)

        results = {
            "rescue": self.rescue_name,
            "drug": self.drug_name,
            "n_conformers": len(conformers),
            "successful_dockings": successful_dockings,
            "affinities": affinities.tolist(),
            "consensus_affinity": float(consensus_affinity),
            "affinity_std": float(affinity_std),
            "best_affinity": float(affinities.min()),
            "worst_affinity": float(affinities.max())
        }

        # Save results
        results_file = self.output_dir / "docking_results.json"
        with open(results_file, 'w') as f:
            json.dump(results, f, indent=2)

        print("\n" + "=" * 80)
        print("✅ Ensemble Docking Complete")
        print("=" * 80)
        print(f"Rescue: {self.rescue_name}")
        print(f"Drug: {self.drug_name}")
        print(f"Successful dockings: {successful_dockings}/{len(conformers)}")
        print(f"\nResults:")
        print(f"  - Consensus affinity: {consensus_affinity:.2f} ± {affinity_std:.2f} kcal/mol")
        print(f"  - Best affinity: {affinities.min():.2f} kcal/mol")
        print(f"  - Worst affinity: {affinities.max():.2f} kcal/mol")
        print(f"\nOutput:")
        print(f"  - Results: {results_file}")
        print(f"  - Docked structures: {self.output_dir}")

        return results


def main():
    parser = argparse.ArgumentParser(description="Run ensemble docking")
    parser.add_argument("--rescue", required=True, help="Rescue name (e.g., 'A189S_M133L_S95T')")
    parser.add_argument("--drug", required=True, help="Drug name (e.g., 'COTI-2')")

    args = parser.parse_args()

    # Create docking instance
    docking = EnsembleDocking(args.rescue, args.drug)

    # Run docking
    try:
        results = docking.run_ensemble_docking()
        exit(0)
    except Exception as e:
        print(f"\n❌ Docking pipeline failed: {e}")
        import traceback
        traceback.print_exc()
        exit(1)


if __name__ == "__main__":
    main()
