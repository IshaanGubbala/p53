"""Validate protein structures for consistency in residue numbering and chains."""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

from Bio import PDB


logger = logging.getLogger(__name__)


def extract_core_domain(
    input_pdb: Path,
    output_pdb: Path,
    chain_id: str,
    start_res: int,
    end_res: int,
) -> bool:
    """Extract core domain from PDB structure.

    Args:
        input_pdb: Input PDB file
        output_pdb: Output PDB file (core domain only)
        chain_id: Chain identifier
        start_res: Start residue number (inclusive)
        end_res: End residue number (inclusive)

    Returns:
        True if successful, False otherwise
    """
    try:
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("protein", input_pdb)

        # Create new structure with only core domain
        io = PDB.PDBIO()

        class CoreDomainSelect(PDB.Select):
            def accept_residue(self, residue):
                if residue.get_parent().id != chain_id:
                    return 0
                res_id = residue.id[1]
                return start_res <= res_id <= end_res

        io.set_structure(structure)
        output_pdb.parent.mkdir(parents=True, exist_ok=True)
        io.save(str(output_pdb), CoreDomainSelect())

        logger.info(f"Extracted core domain {start_res}-{end_res} (chain {chain_id}) from {input_pdb.name}")
        return True

    except Exception as e:
        logger.error(f"Failed to extract core domain: {e}")
        return False


def validate_residue_numbering(
    pdb_path: Path,
    chain_id: str,
    expected_start: int,
    expected_end: int,
) -> dict[str, any]:
    """Validate residue numbering in PDB structure.

    Args:
        pdb_path: Path to PDB file
        chain_id: Chain identifier
        expected_start: Expected first residue number
        expected_end: Expected last residue number

    Returns:
        Dictionary with validation results
    """
    try:
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("protein", pdb_path)

        # Get chain
        model = structure[0]
        if chain_id not in model:
            return {
                "valid": False,
                "reason": f"Chain {chain_id} not found",
                "available_chains": [c.id for c in model.get_chains()],
            }

        chain = model[chain_id]

        # Get residue numbers
        residues = [res for res in chain.get_residues() if PDB.is_aa(res)]

        if not residues:
            return {
                "valid": False,
                "reason": "No amino acid residues found",
            }

        res_numbers = sorted([res.id[1] for res in residues])
        actual_start = res_numbers[0]
        actual_end = res_numbers[-1]

        # Check numbering
        numbering_matches = (actual_start == expected_start and actual_end == expected_end)

        # Check for gaps
        expected_range = set(range(expected_start, expected_end + 1))
        actual_range = set(res_numbers)
        missing_residues = expected_range - actual_range
        extra_residues = actual_range - expected_range

        return {
            "valid": numbering_matches and len(missing_residues) == 0,
            "actual_start": actual_start,
            "actual_end": actual_end,
            "expected_start": expected_start,
            "expected_end": expected_end,
            "n_residues": len(residues),
            "missing_residues": sorted(missing_residues),
            "extra_residues": sorted(extra_residues),
            "numbering_matches": numbering_matches,
        }

    except Exception as e:
        return {
            "valid": False,
            "reason": f"Validation failed: {str(e)}",
        }


def compare_structures(
    pdb1: Path,
    pdb2: Path,
    chain1: str,
    chain2: str,
) -> dict[str, any]:
    """Compare two structures for compatibility.

    Args:
        pdb1: First PDB file
        pdb2: Second PDB file
        chain1: Chain ID in first structure
        chain2: Chain ID in second structure

    Returns:
        Dictionary with comparison results
    """
    try:
        parser = PDB.PDBParser(QUIET=True)

        # Parse structures
        struct1 = parser.get_structure("s1", pdb1)
        struct2 = parser.get_structure("s2", pdb2)

        # Get chains
        chain_obj1 = struct1[0][chain1]
        chain_obj2 = struct2[0][chain2]

        # Get residue numbers
        res1 = [res.id[1] for res in chain_obj1.get_residues() if PDB.is_aa(res)]
        res2 = [res.id[1] for res in chain_obj2.get_residues() if PDB.is_aa(res)]

        # Find common residues
        common_res = sorted(set(res1) & set(res2))

        # Check sequence identity at common positions
        parser1_seq = []
        parser2_seq = []

        for res_num in common_res:
            # Get 3-letter code
            try:
                aa1 = chain_obj1[res_num].resname
                aa2 = chain_obj2[res_num].resname
                parser1_seq.append(aa1)
                parser2_seq.append(aa2)
            except KeyError:
                continue

        # Count matches
        matches = sum(1 for a1, a2 in zip(parser1_seq, parser2_seq) if a1 == a2)
        identity = matches / len(common_res) if common_res else 0.0

        return {
            "compatible": identity > 0.95,  # 95% sequence identity threshold
            "n_common_residues": len(common_res),
            "sequence_identity": identity,
            "pdb1_range": (min(res1), max(res1)),
            "pdb2_range": (min(res2), max(res2)),
            "only_in_pdb1": sorted(set(res1) - set(res2)),
            "only_in_pdb2": sorted(set(res2) - set(res1)),
        }

    except Exception as e:
        return {
            "compatible": False,
            "reason": f"Comparison failed: {str(e)}",
        }


def get_structure_info(pdb_path: Path, chain_id: str) -> dict[str, any]:
    """Get basic structure information.

    Args:
        pdb_path: Path to PDB file
        chain_id: Chain identifier

    Returns:
        Dictionary with structure info
    """
    try:
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("protein", pdb_path)

        model = structure[0]
        chain = model[chain_id]

        residues = [res for res in chain.get_residues() if PDB.is_aa(res)]
        res_numbers = sorted([res.id[1] for res in residues])

        return {
            "pdb_path": str(pdb_path),
            "chain_id": chain_id,
            "n_residues": len(residues),
            "start_residue": res_numbers[0] if res_numbers else None,
            "end_residue": res_numbers[-1] if res_numbers else None,
            "residue_range": (min(res_numbers), max(res_numbers)) if res_numbers else None,
        }

    except Exception as e:
        return {
            "pdb_path": str(pdb_path),
            "error": str(e),
        }
