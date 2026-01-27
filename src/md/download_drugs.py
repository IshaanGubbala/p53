#!/usr/bin/env python3
"""
Download and prepare drug structures from PubChem.

Usage:
    python src/md/download_drugs.py

Requires:
    - pubchempy
    - rdkit
    - openbabel (system)
"""

import json
import subprocess
from pathlib import Path
from typing import Dict, Any

try:
    import pubchempy as pcp
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    print("❌ Missing dependencies. Install with:")
    print("   pip install pubchempy rdkit")
    exit(1)


def load_drug_metadata() -> Dict[str, Any]:
    """Load drug metadata from JSON file."""
    metadata_path = Path("Data/raw/drugs/metadata.json")
    with open(metadata_path) as f:
        data = json.load(f)
    return data["drugs"]


def download_from_pubchem(drug_name: str, cid: int, output_dir: Path) -> bool:
    """
    Download 3D structure from PubChem.

    Args:
        drug_name: Name of the drug
        cid: PubChem compound ID
        output_dir: Directory to save structure

    Returns:
        True if successful, False otherwise
    """
    print(f"\n📥 Downloading {drug_name} (CID: {cid}) from PubChem...")

    try:
        # Get compound
        compounds = pcp.get_compounds(cid, 'cid')
        if not compounds:
            print(f"  ❌ No compound found for CID {cid}")
            return False

        compound = compounds[0]

        # Get 3D conformer (SDF format)
        sdf_path = output_dir / f"{drug_name}.sdf"
        pcp.download('SDF', sdf_path, cid, 'cid', overwrite=True, record_type='3d')

        if not sdf_path.exists():
            print(f"  ❌ Failed to download SDF for {drug_name}")
            return False

        print(f"  ✅ Downloaded {drug_name}.sdf")
        return True

    except Exception as e:
        print(f"  ❌ Error downloading {drug_name}: {e}")
        return False


def generate_3d_from_smiles(drug_name: str, smiles: str, output_dir: Path) -> bool:
    """
    Generate 3D structure from SMILES using RDKit.

    Args:
        drug_name: Name of the drug
        smiles: SMILES string
        output_dir: Directory to save structure

    Returns:
        True if successful, False otherwise
    """
    print(f"\n🧪 Generating 3D structure for {drug_name} from SMILES...")

    try:
        # Parse SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"  ❌ Invalid SMILES: {smiles}")
            return False

        # Add hydrogens
        mol = Chem.AddHs(mol)

        # Generate 3D coordinates
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol, maxIters=500)

        # Save to SDF
        sdf_path = output_dir / f"{drug_name}.sdf"
        writer = Chem.SDWriter(str(sdf_path))
        writer.write(mol)
        writer.close()

        print(f"  ✅ Generated {drug_name}.sdf from SMILES")
        return True

    except Exception as e:
        print(f"  ❌ Error generating 3D structure: {e}")
        return False


def convert_to_pdb(drug_name: str, input_dir: Path, output_dir: Path) -> bool:
    """
    Convert SDF to PDB using Open Babel.

    Args:
        drug_name: Name of the drug
        input_dir: Directory containing SDF file
        output_dir: Directory to save PDB

    Returns:
        True if successful, False otherwise
    """
    sdf_path = input_dir / f"{drug_name}.sdf"
    pdb_path = output_dir / f"{drug_name}.pdb"

    if not sdf_path.exists():
        print(f"  ❌ SDF file not found: {sdf_path}")
        return False

    print(f"  🔄 Converting {drug_name}.sdf to PDB...")

    try:
        # Use Open Babel for conversion
        cmd = [
            "obabel",
            str(sdf_path),
            "-O", str(pdb_path),
            "-h",  # Add hydrogens
            "--gen3d"  # Generate 3D coordinates if needed
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, check=True)

        if pdb_path.exists():
            print(f"  ✅ Converted to {drug_name}.pdb")
            return True
        else:
            print(f"  ❌ Conversion failed")
            return False

    except subprocess.CalledProcessError as e:
        print(f"  ⚠️  Open Babel not available, skipping conversion")
        print(f"     (Install with: conda install -c conda-forge openbabel)")
        return False
    except Exception as e:
        print(f"  ❌ Error during conversion: {e}")
        return False


def minimize_structure(drug_name: str, pdb_dir: Path) -> bool:
    """
    Energy minimize structure using Open Babel.

    Args:
        drug_name: Name of the drug
        pdb_dir: Directory containing PDB file

    Returns:
        True if successful, False otherwise
    """
    pdb_path = pdb_dir / f"{drug_name}.pdb"
    minimized_path = pdb_dir / f"{drug_name}_minimized.pdb"

    if not pdb_path.exists():
        return False

    print(f"  ⚡ Energy minimizing {drug_name}...")

    try:
        cmd = [
            "obminimize",
            "-ff", "MMFF94",
            "-n", "500",
            str(pdb_path)
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, check=True)

        # Save minimized structure
        with open(minimized_path, 'w') as f:
            f.write(result.stdout)

        print(f"  ✅ Minimized to {drug_name}_minimized.pdb")
        return True

    except subprocess.CalledProcessError:
        print(f"  ⚠️  Energy minimization skipped (obminimize not available)")
        return False
    except Exception as e:
        print(f"  ❌ Error during minimization: {e}")
        return False


def prepare_for_docking(drug_name: str, pdb_dir: Path, output_dir: Path) -> bool:
    """
    Prepare drug for docking (convert to PDBQT format).

    Args:
        drug_name: Name of the drug
        pdb_dir: Directory containing PDB file
        output_dir: Directory to save PDBQT

    Returns:
        True if successful, False otherwise
    """
    # Use minimized if available, otherwise original
    minimized_path = pdb_dir / f"{drug_name}_minimized.pdb"
    pdb_path = minimized_path if minimized_path.exists() else pdb_dir / f"{drug_name}.pdb"
    pdbqt_path = output_dir / f"{drug_name}.pdbqt"

    if not pdb_path.exists():
        return False

    print(f"  🎯 Preparing {drug_name} for docking...")

    try:
        # Use AutoDockTools prepare_ligand4.py
        cmd = [
            "prepare_ligand4.py",
            "-l", str(pdb_path),
            "-o", str(pdbqt_path),
            "-A", "hydrogens"  # Add hydrogens
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, check=True)

        if pdbqt_path.exists():
            print(f"  ✅ Prepared {drug_name}.pdbqt for docking")
            return True
        else:
            print(f"  ❌ PDBQT conversion failed")
            return False

    except subprocess.CalledProcessError:
        print(f"  ⚠️  PDBQT conversion skipped (prepare_ligand4.py not available)")
        print(f"     (Install AutoDockTools from MGLTools)")
        return False
    except Exception as e:
        print(f"  ❌ Error preparing for docking: {e}")
        return False


def main():
    """Main workflow: download and prepare all drugs."""

    print("=" * 80)
    print("Drug Structure Download & Preparation Pipeline")
    print("=" * 80)

    # Setup directories
    base_dir = Path("Data/raw/drugs")
    structures_dir = base_dir / "structures"
    prepared_dir = base_dir / "prepared"

    structures_dir.mkdir(parents=True, exist_ok=True)
    prepared_dir.mkdir(parents=True, exist_ok=True)

    # Load drug metadata
    drugs = load_drug_metadata()

    print(f"\n📋 Found {len(drugs)} drugs to process\n")

    # Process each drug
    results = {"success": [], "failed": []}

    for drug_name, drug_info in drugs.items():
        print(f"\n{'─' * 80}")
        print(f"Processing: {drug_name}")
        print(f"{'─' * 80}")

        success = False

        # Try PubChem download first
        cid = drug_info.get("pubchem_cid")
        if cid:
            success = download_from_pubchem(drug_name, cid, structures_dir)

        # Fallback to SMILES if PubChem fails or no CID
        if not success and "smiles" in drug_info:
            success = generate_3d_from_smiles(
                drug_name,
                drug_info["smiles"],
                structures_dir
            )

        if not success:
            print(f"  ❌ Failed to obtain structure for {drug_name}")
            results["failed"].append(drug_name)
            continue

        # Convert to PDB
        convert_to_pdb(drug_name, structures_dir, structures_dir)

        # Energy minimize
        minimize_structure(drug_name, structures_dir)

        # Prepare for docking
        prepare_for_docking(drug_name, structures_dir, prepared_dir)

        results["success"].append(drug_name)

    # Summary
    print("\n" + "=" * 80)
    print("Summary")
    print("=" * 80)
    print(f"✅ Successfully processed: {len(results['success'])} drugs")
    for drug in results["success"]:
        print(f"   - {drug}")

    if results["failed"]:
        print(f"\n❌ Failed to process: {len(results['failed'])} drugs")
        for drug in results["failed"]:
            print(f"   - {drug}")

    print("\n📁 Output directories:")
    print(f"   - Structures: {structures_dir}")
    print(f"   - Prepared (PDBQT): {prepared_dir}")
    print()


if __name__ == "__main__":
    main()
