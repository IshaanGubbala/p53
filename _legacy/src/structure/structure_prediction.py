"""
Structure Prediction Module for p53 Rescue Mutants

IMPORTANT: This module provides proper structure prediction for mutant sequences.
DO NOT use EvoEF2 BuildMutant for structure generation - it only swaps side chains
without backbone relaxation, creating strained geometries that cause 6-8 Å RMSD
drift in MD simulations.

Use this module for:
- Generating mutant structures for MD simulation
- Accurate structural analysis

Use EvoEF2 for:
- ΔΔG scoring (stability prediction)
- Quick mutation screening
- Ranking rescue candidates

Supported methods:
1. ESMFold (fast, local) - recommended for quick predictions
2. ColabFold/AlphaFold2 (accurate) - recommended for final structures
3. MODELLER (homology) - for multi-domain modeling
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np

logger = logging.getLogger(__name__)

# p53 wildtype sequence (canonical, UniProt P04637)
P53_WILDTYPE_SEQUENCE = (
    "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP"
    "DEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAK"
    "SVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHE"
    "RCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNS"
    "SCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELP"
    "PGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPG"
    "GSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD"
)


@dataclass
class StructureValidation:
    """Results of structure validation."""
    is_valid: bool
    rmsd_from_template: Optional[float]  # Å
    clash_score: float  # Lower is better
    ramachandran_favored: float  # Percentage
    ramachandran_outliers: float  # Percentage
    warnings: list[str]
    errors: list[str]


@dataclass
class PredictionResult:
    """Result of structure prediction."""
    pdb_path: Path
    method: str
    confidence: Optional[float]  # pLDDT for AlphaFold/ESMFold
    validation: Optional[StructureValidation]
    mutations: list[str]
    sequence: str


def apply_mutations_to_sequence(
    sequence: str,
    mutations: list[str],
) -> str:
    """
    Apply mutations to a sequence.

    Args:
        sequence: Wildtype sequence
        mutations: List of mutations in format "X123Y" (1-indexed)

    Returns:
        Mutated sequence

    Raises:
        ValueError: If mutation is invalid or doesn't match sequence
    """
    import re

    seq_list = list(sequence)

    for mut in mutations:
        match = re.match(r'^([A-Z])(\d+)([A-Z])$', mut.upper())
        if not match:
            raise ValueError(f"Invalid mutation format: {mut}. Expected 'X123Y'")

        ref, pos_str, alt = match.groups()
        pos = int(pos_str) - 1  # Convert to 0-indexed

        if pos < 0 or pos >= len(seq_list):
            raise ValueError(f"Position {pos_str} out of range (1-{len(seq_list)})")

        if seq_list[pos] != ref:
            raise ValueError(
                f"Reference mismatch at position {pos_str}: "
                f"expected {ref}, found {seq_list[pos]}"
            )

        seq_list[pos] = alt
        logger.debug(f"Applied mutation {mut}: {ref}{pos_str} -> {alt}")

    return ''.join(seq_list)


def predict_structure_esmfold(
    sequence: str,
    output_path: Path,
    device: str = "cuda",
) -> PredictionResult:
    """
    Predict structure using ESMFold (fast, local prediction).

    ESMFold is a single-sequence structure prediction model that runs locally
    without requiring MSA computation. It's faster than AlphaFold2 but slightly
    less accurate.

    Args:
        sequence: Amino acid sequence
        output_path: Path to save predicted PDB
        device: "cuda" or "cpu"

    Returns:
        PredictionResult with predicted structure

    Note:
        Requires: pip install fair-esm
        GPU with 16GB+ VRAM recommended for full p53 (393 aa)
    """
    try:
        import esm
        import torch
    except ImportError:
        raise ImportError(
            "ESMFold requires fair-esm package. Install with: pip install fair-esm"
        )

    logger.info(f"Loading ESMFold model on {device}...")
    model = esm.pretrained.esmfold_v1()
    model = model.eval()

    if device == "cuda" and torch.cuda.is_available():
        model = model.cuda()
        logger.info(f"Using GPU: {torch.cuda.get_device_name()}")
    else:
        device = "cpu"
        logger.info("Using CPU (slow)")

    # Predict structure
    logger.info(f"Predicting structure for {len(sequence)} residues...")
    with torch.no_grad():
        output = model.infer_pdb(sequence)

    # Extract pLDDT confidence
    # ESMFold returns pLDDT in B-factor column
    plddt_values = []
    for line in output.split('\n'):
        if line.startswith('ATOM'):
            try:
                bfactor = float(line[60:66].strip())
                plddt_values.append(bfactor)
            except (ValueError, IndexError):
                pass

    mean_plddt = np.mean(plddt_values) if plddt_values else None

    # Save PDB
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(output)
    logger.info(f"Saved structure to {output_path}")

    if mean_plddt:
        logger.info(f"Mean pLDDT confidence: {mean_plddt:.1f}")

    return PredictionResult(
        pdb_path=output_path,
        method="ESMFold",
        confidence=mean_plddt,
        validation=None,
        mutations=[],
        sequence=sequence,
    )


def predict_structure_colabfold(
    sequence: str,
    output_dir: Path,
    name: str = "mutant",
    num_recycle: int = 3,
    use_amber_relax: bool = True,
) -> PredictionResult:
    """
    Predict structure using ColabFold (AlphaFold2 with MMseqs2 MSA).

    This is the most accurate method but requires internet access for MSA
    computation via ColabFold servers.

    Args:
        sequence: Amino acid sequence
        output_dir: Directory for output files
        name: Name prefix for output files
        num_recycle: Number of recycling iterations (3-5 recommended)
        use_amber_relax: Whether to run AMBER relaxation

    Returns:
        PredictionResult with predicted structure

    Note:
        Requires: pip install colabfold
        Or run in Google Colab notebook
    """
    try:
        from colabfold.batch import run as colabfold_run
    except ImportError:
        raise ImportError(
            "ColabFold not installed. Install with: "
            "pip install 'colabfold[alphafold] @ git+https://github.com/sokrypton/ColabFold'"
        )

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Write FASTA
    fasta_path = output_dir / f"{name}.fasta"
    fasta_path.write_text(f">{name}\n{sequence}\n")

    logger.info(f"Running ColabFold for {len(sequence)} residues...")

    # Run ColabFold
    colabfold_run(
        queries=str(fasta_path),
        results_dir=str(output_dir),
        num_recycle=num_recycle,
        use_amber=use_amber_relax,
        num_models=1,
    )

    # Find best ranked model
    relaxed_models = list(output_dir.glob("*_relaxed_rank_001*.pdb"))
    unrelaxed_models = list(output_dir.glob("*_unrelaxed_rank_001*.pdb"))

    if relaxed_models:
        pdb_path = relaxed_models[0]
    elif unrelaxed_models:
        pdb_path = unrelaxed_models[0]
    else:
        raise RuntimeError(f"No output PDB found in {output_dir}")

    # Read pLDDT from JSON
    json_files = list(output_dir.glob("*_scores_rank_001*.json"))
    mean_plddt = None
    if json_files:
        with open(json_files[0]) as f:
            scores = json.load(f)
            if 'plddt' in scores:
                mean_plddt = np.mean(scores['plddt'])

    logger.info(f"ColabFold prediction complete: {pdb_path}")
    if mean_plddt:
        logger.info(f"Mean pLDDT: {mean_plddt:.1f}")

    return PredictionResult(
        pdb_path=pdb_path,
        method="ColabFold",
        confidence=mean_plddt,
        validation=None,
        mutations=[],
        sequence=sequence,
    )


def validate_structure(
    pdb_path: Path,
    template_pdb: Optional[Path] = None,
    check_clashes: bool = True,
    check_ramachandran: bool = True,
) -> StructureValidation:
    """
    Validate a predicted structure.

    Args:
        pdb_path: Path to PDB file to validate
        template_pdb: Optional template for RMSD calculation
        check_clashes: Whether to check for steric clashes
        check_ramachandran: Whether to check Ramachandran angles

    Returns:
        StructureValidation with results
    """
    warnings = []
    errors = []

    # Read PDB
    try:
        from Bio.PDB import PDBParser, PDBIO
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("mutant", str(pdb_path))
    except ImportError:
        warnings.append("BioPython not installed - limited validation")
        return StructureValidation(
            is_valid=True,
            rmsd_from_template=None,
            clash_score=0.0,
            ramachandran_favored=0.0,
            ramachandran_outliers=0.0,
            warnings=warnings,
            errors=errors,
        )

    # RMSD from template
    rmsd = None
    if template_pdb and template_pdb.exists():
        try:
            from Bio.PDB import Superimposer
            template = parser.get_structure("template", str(template_pdb))

            # Get CA atoms
            mut_ca = [atom for atom in structure.get_atoms() if atom.name == "CA"]
            tmpl_ca = [atom for atom in template.get_atoms() if atom.name == "CA"]

            if len(mut_ca) == len(tmpl_ca):
                sup = Superimposer()
                sup.set_atoms(tmpl_ca, mut_ca)
                rmsd = sup.rms

                if rmsd > 5.0:
                    warnings.append(f"High RMSD from template: {rmsd:.2f} Å")
                elif rmsd > 3.0:
                    warnings.append(f"Moderate RMSD from template: {rmsd:.2f} Å")
            else:
                warnings.append(
                    f"Chain length mismatch: mutant={len(mut_ca)}, template={len(tmpl_ca)}"
                )
        except Exception as e:
            warnings.append(f"RMSD calculation failed: {e}")

    # Clash detection
    clash_score = 0.0
    if check_clashes:
        try:
            atoms = list(structure.get_atoms())
            n_atoms = len(atoms)
            n_clashes = 0

            for i in range(n_atoms):
                for j in range(i + 2, n_atoms):  # Skip bonded atoms
                    dist = atoms[i] - atoms[j]
                    if dist < 1.5:  # Very close contact
                        n_clashes += 1

            clash_score = n_clashes / max(n_atoms, 1)

            if clash_score > 0.01:
                warnings.append(f"High clash score: {clash_score:.4f} ({n_clashes} clashes)")
        except Exception as e:
            warnings.append(f"Clash detection failed: {e}")

    # Ramachandran (simplified)
    rama_favored = 0.0
    rama_outliers = 0.0
    if check_ramachandran:
        try:
            from Bio.PDB.Polypeptide import PPBuilder
            ppb = PPBuilder()

            n_residues = 0
            n_favored = 0
            n_outliers = 0

            for pp in ppb.build_peptides(structure):
                phi_psi = pp.get_phi_psi_list()
                for phi, psi in phi_psi:
                    if phi is None or psi is None:
                        continue
                    n_residues += 1

                    # Simplified Ramachandran check
                    # Favored regions (very simplified)
                    phi_deg = np.degrees(phi)
                    psi_deg = np.degrees(psi)

                    # Alpha helix: phi ~ -60, psi ~ -45
                    # Beta sheet: phi ~ -120, psi ~ 120
                    in_favored = (
                        (-180 < phi_deg < -30 and -90 < psi_deg < 45) or  # Alpha
                        (-180 < phi_deg < -60 and 90 < psi_deg < 180) or  # Beta
                        (-180 < phi_deg < -60 and -180 < psi_deg < -120)  # Beta
                    )

                    if in_favored:
                        n_favored += 1
                    elif abs(phi_deg) > 160 or abs(psi_deg) > 160:
                        n_outliers += 1

            if n_residues > 0:
                rama_favored = 100 * n_favored / n_residues
                rama_outliers = 100 * n_outliers / n_residues

            if rama_outliers > 5.0:
                warnings.append(f"High Ramachandran outliers: {rama_outliers:.1f}%")
        except Exception as e:
            warnings.append(f"Ramachandran check failed: {e}")

    # Overall validity
    is_valid = len(errors) == 0 and clash_score < 0.05

    return StructureValidation(
        is_valid=is_valid,
        rmsd_from_template=rmsd,
        clash_score=clash_score,
        ramachandran_favored=rama_favored,
        ramachandran_outliers=rama_outliers,
        warnings=warnings,
        errors=errors,
    )


def predict_mutant_structure(
    mutations: list[str],
    output_path: Path,
    method: str = "esmfold",
    template_pdb: Optional[Path] = None,
    validate: bool = True,
    cache_dir: Optional[Path] = None,
    wildtype_sequence: str = P53_WILDTYPE_SEQUENCE,
) -> PredictionResult:
    """
    Predict structure of a p53 mutant.

    This is the main entry point for structure prediction. Use this instead of
    EvoEF2 BuildMutant for generating structures for MD simulation or detailed
    structural analysis.

    Args:
        mutations: List of mutations in "X123Y" format (1-indexed)
        output_path: Where to save the predicted PDB
        method: "esmfold" (fast, local) or "colabfold" (accurate, requires internet)
        template_pdb: Optional template for validation RMSD
        validate: Whether to run structure validation
        cache_dir: Optional cache directory for reuse
        wildtype_sequence: Base sequence (default: p53)

    Returns:
        PredictionResult with predicted structure and validation

    Example:
        >>> result = predict_mutant_structure(
        ...     mutations=["S95T", "M133L", "A189S"],
        ...     output_path=Path("mutant.pdb"),
        ...     method="esmfold"
        ... )
        >>> print(f"Saved to {result.pdb_path}, pLDDT={result.confidence:.1f}")
    """
    # Check cache
    if cache_dir:
        cache_key = hashlib.md5(
            f"{','.join(sorted(mutations))}_{method}".encode()
        ).hexdigest()[:12]
        cached_pdb = cache_dir / f"structure_{cache_key}.pdb"

        if cached_pdb.exists():
            logger.info(f"Using cached structure: {cached_pdb}")
            shutil.copy(cached_pdb, output_path)

            validation = None
            if validate:
                validation = validate_structure(output_path, template_pdb)

            return PredictionResult(
                pdb_path=output_path,
                method=f"{method} (cached)",
                confidence=None,
                validation=validation,
                mutations=mutations,
                sequence=apply_mutations_to_sequence(wildtype_sequence, mutations),
            )

    # Apply mutations to sequence
    mutant_sequence = apply_mutations_to_sequence(wildtype_sequence, mutations)
    logger.info(f"Predicting structure for {len(mutations)} mutation(s): {', '.join(mutations)}")

    # Predict structure
    if method.lower() == "esmfold":
        result = predict_structure_esmfold(mutant_sequence, output_path)
    elif method.lower() in ("colabfold", "alphafold"):
        result = predict_structure_colabfold(
            mutant_sequence,
            output_path.parent,
            name=output_path.stem,
        )
        # Copy to desired location if different
        if result.pdb_path != output_path:
            shutil.copy(result.pdb_path, output_path)
            result.pdb_path = output_path
    else:
        raise ValueError(f"Unknown method: {method}. Use 'esmfold' or 'colabfold'")

    result.mutations = mutations
    result.sequence = mutant_sequence

    # Validate
    if validate:
        result.validation = validate_structure(output_path, template_pdb)

        if result.validation.warnings:
            for warn in result.validation.warnings:
                logger.warning(warn)
        if result.validation.errors:
            for err in result.validation.errors:
                logger.error(err)

    # Cache result
    if cache_dir:
        cache_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy(output_path, cached_pdb)
        logger.debug(f"Cached structure to {cached_pdb}")

    return result


# Convenience function for MD simulation
def prepare_structure_for_md(
    mutations: list[str],
    output_dir: Path,
    method: str = "esmfold",
) -> Path:
    """
    Prepare a mutant structure for MD simulation.

    This function predicts the structure and validates it's suitable for MD.

    Args:
        mutations: List of mutations
        output_dir: Output directory
        method: Prediction method

    Returns:
        Path to validated PDB file ready for MD

    Raises:
        RuntimeError: If structure validation fails
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    mutation_str = "_".join(mutations) if mutations else "wildtype"
    pdb_path = output_dir / f"p53_{mutation_str}.pdb"

    result = predict_mutant_structure(
        mutations=mutations,
        output_path=pdb_path,
        method=method,
        validate=True,
    )

    if result.validation and not result.validation.is_valid:
        raise RuntimeError(
            f"Structure validation failed: {result.validation.errors}"
        )

    if result.confidence and result.confidence < 50:
        logger.warning(
            f"Low prediction confidence ({result.confidence:.1f}). "
            "Consider using ColabFold for better accuracy."
        )

    return pdb_path
