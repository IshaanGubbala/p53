"""
Tetramer interface analysis for rescue mutations.

This module identifies rescue mutations that might disrupt p53 tetramerization,
causing dominant-negative effects.

CRITICAL: p53 functions as a tetramer. A mutation that stabilizes the monomer
but disrupts tetramer formation can "poison" the entire complex.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, NeighborSearch

from src.core.logging import get_logger

logger = get_logger(__name__)


# Core domain tetramer interface residues (from 3KMD structure analysis)
CORE_DOMAIN_INTERFACE_RESIDUES = [
    93, 94, 95, 96, 98, 99, 100, 101,  # N-terminal interface
    140,  # Loop region
    166, 167, 170,  # DNA-binding region
    176, 177, 178, 179, 180, 181,  # Zinc-binding region
    198, 199, 200, 201, 202,  # H2 helix
    212,  # Loop
    224, 225, 226, 231, 233,  # L3 loop, S7-S8
    242, 243, 244,  # Zinc-binding region
    267,  # H2 helix
]


def identify_tetramer_interface(
    pdb_path: Path,
    protein_chains: list[str] = ["A", "B", "C", "D"],
    distance_threshold: float = 5.0,
) -> dict[int, dict[str, Any]]:
    """
    Identify residues at tetramer interfaces.

    Args:
        pdb_path: Path to PDB structure with p53 tetramer
        protein_chains: List of protein chain IDs
        distance_threshold: Maximum distance (Å) for interface contact

    Returns:
        Dictionary mapping residue position to interface information
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("tetramer", str(pdb_path))
    model = structure[0]

    # Get atoms per chain
    chain_atoms = {}
    for chain in model:
        chain_id = chain.get_id()
        if chain_id in protein_chains:
            atoms = []
            for residue in chain:
                if residue.get_id()[0] == " ":  # Standard residue
                    atoms.extend(residue.get_atoms())
            chain_atoms[chain_id] = atoms

    logger.info(f"Found {len(chain_atoms)} protein chains")

    # Find interface residues
    interface_contacts = {}

    for chain1_id, atoms1 in chain_atoms.items():
        for chain2_id, atoms2 in chain_atoms.items():
            if chain2_id <= chain1_id:  # Avoid duplicates
                continue

            # Build neighbor search for chain2
            chain2_search = NeighborSearch(atoms2)

            # Find chain1 atoms near chain2
            for atom1 in atoms1:
                nearby = chain2_search.search(atom1.coord, distance_threshold)
                if nearby:
                    res_id = atom1.get_parent().get_id()[1]
                    res_name = atom1.get_parent().get_resname()

                    if res_id not in interface_contacts:
                        interface_contacts[res_id] = {
                            "position": int(res_id),
                            "residue_name": str(res_name),
                            "n_contacts": 0,
                            "min_distance": float('inf'),
                            "contact_pairs": [],
                        }

                    # Calculate distance
                    for atom2 in nearby:
                        distance = float(np.linalg.norm(atom1.coord - atom2.coord))
                        if distance <= distance_threshold:
                            interface_contacts[res_id]["n_contacts"] += 1
                            interface_contacts[res_id]["min_distance"] = float(min(
                                interface_contacts[res_id]["min_distance"], distance
                            ))
                            interface_contacts[res_id]["contact_pairs"].append({
                                "chain1": chain1_id,
                                "chain2": chain2_id,
                                "distance": distance,
                            })

    logger.info(f"Identified {len(interface_contacts)} interface residues")

    return interface_contacts


def check_tetramer_interface_conflict(
    positions: list[int],
    proximity_threshold: float = 10.0,
) -> dict[str, Any]:
    """
    Check if rescue positions overlap with or are near tetramer interfaces.

    Args:
        positions: List of mutation positions
        proximity_threshold: Distance (Å) to consider "near" interface (default: 10.0)

    Returns:
        Dict with conflict assessment
    """
    interface_positions = set(CORE_DOMAIN_INTERFACE_RESIDUES)

    direct_conflicts = [pos for pos in positions if pos in interface_positions]

    # For proximity, we'd need structure coordinates - for now just flag direct hits
    warnings = []
    risk_level = "safe"

    if direct_conflicts:
        warnings.append(f"Direct interface conflict at position(s): {', '.join(map(str, direct_conflicts))}")
        risk_level = "high"

    return {
        "has_interface_conflict": len(direct_conflicts) > 0,
        "interface_positions": direct_conflicts,
        "risk_level": risk_level,
        "warnings": "; ".join(warnings) if warnings else None,
        "n_interface_conflicts": len(direct_conflicts),
    }


def estimate_tetramer_impact(
    mutation: str,
    interface_residues: set[int] | None = None,
) -> dict[str, Any]:
    """
    Estimate the impact of a mutation on tetramer stability.

    This is a heuristic approximation. For rigorous analysis,
    use FoldX InterfaceAnalyzer or Rosetta.

    Args:
        mutation: Mutation string (e.g., "R175H")
        interface_residues: Set of interface positions

    Returns:
        Dict with impact assessment
    """
    if interface_residues is None:
        interface_residues = set(CORE_DOMAIN_INTERFACE_RESIDUES)

    # Parse mutation
    if len(mutation) < 3:
        return {"error": "Invalid mutation format"}

    wt_aa = mutation[0]
    position = int(''.join(filter(str.isdigit, mutation)))
    mut_aa = mutation[-1]

    # Check if position is at interface
    if position not in interface_residues:
        return {
            "position": int(position),
            "at_interface": False,
            "impact": "none",
            "risk_level": "safe",
            "explanation": f"Position {position} is not at tetramer interface",
        }

    # Amino acid properties
    charged_positive = ["R", "K", "H"]
    charged_negative = ["D", "E"]
    polar = ["S", "T", "N", "Q"]
    hydrophobic = ["A", "V", "L", "I", "M", "F", "W", "Y"]
    small = ["G", "A", "S"]

    # Assess mutation impact
    impact_score = 0
    explanations = []

    # 1. Loss of electrostatic interaction
    if wt_aa in charged_positive and mut_aa not in charged_positive:
        impact_score += 2
        explanations.append(f"Loss of positive charge at interface ({wt_aa}→{mut_aa})")
    elif wt_aa in charged_negative and mut_aa not in charged_negative:
        impact_score += 2
        explanations.append(f"Loss of negative charge at interface ({wt_aa}→{mut_aa})")

    # 2. Introduction of steric clash
    if wt_aa in small and mut_aa not in small:
        impact_score += 1.5
        explanations.append(f"Introduction of bulky residue at interface ({wt_aa}→{mut_aa})")

    # 3. Loss of H-bond capability
    if wt_aa in (polar + charged_positive + charged_negative) and mut_aa in hydrophobic:
        impact_score += 1.5
        explanations.append(f"Loss of H-bond capability at interface ({wt_aa}→{mut_aa})")

    # 4. Proline at interface (helix breaker)
    if mut_aa == "P":
        impact_score += 3
        explanations.append("Introduction of proline (helix breaker) at interface")

    # Determine risk level
    if impact_score >= 4:
        risk_level = "critical"
        impact = "severe"
    elif impact_score >= 2.5:
        risk_level = "high"
        impact = "moderate"
    elif impact_score >= 1:
        risk_level = "medium"
        impact = "mild"
    else:
        risk_level = "low"
        impact = "minimal"

    return {
        "position": int(position),
        "mutation": str(mutation),
        "wt_aa": str(wt_aa),
        "mut_aa": str(mut_aa),
        "at_interface": True,
        "impact_score": float(round(impact_score, 2)),
        "impact": str(impact),
        "risk_level": str(risk_level),
        "explanations": str("; ".join(explanations) if explanations else "Minimal impact expected"),
    }


def analyze_rescue_tetramer_impact(
    rescues_df: pd.DataFrame,
    tetramer_pdb: Path | None = None,
    filter_high_risk: bool = True,
) -> pd.DataFrame:
    """
    Analyze tetramer interface impact for all rescue mutations.

    Args:
        rescues_df: DataFrame with rescue candidates
        tetramer_pdb: Optional path to tetramer structure for detailed analysis
        filter_high_risk: If True, flag rescues with high tetramer risk

    Returns:
        DataFrame with tetramer impact columns
    """
    interface_residues = set(CORE_DOMAIN_INTERFACE_RESIDUES)

    if tetramer_pdb and tetramer_pdb.exists():
        logger.info("Loading tetramer structure for detailed interface analysis...")
        interface_contacts = identify_tetramer_interface(tetramer_pdb)
        interface_residues = set(interface_contacts.keys())

    logger.info(f"Using {len(interface_residues)} interface residues")

    # Analyze each rescue
    tetramer_results = []

    for idx, row in rescues_df.iterrows():
        rescue_mutations = row.get("rescue_mutations", "")

        if pd.isna(rescue_mutations):
            continue

        # Parse mutations
        mutations = [m.strip() for m in str(rescue_mutations).split(",")]

        # Check each mutation for tetramer interface impact
        max_impact_score = 0
        max_risk_level = "safe"
        impacts = []

        for mut in mutations:
            if len(mut) < 3:
                continue

            impact = estimate_tetramer_impact(mut, interface_residues)
            impacts.append(impact)

            if impact.get("impact_score", 0) > max_impact_score:
                max_impact_score = impact["impact_score"]
                max_risk_level = impact["risk_level"]

        tetramer_results.append({
            "rescue_mutations": rescue_mutations,
            "tetramer_impact_score": max_impact_score,
            "tetramer_risk": max_risk_level,
            "at_interface_any": any(imp.get("at_interface", False) for imp in impacts),
            "tetramer_details": impacts,
        })

    # Merge with original DataFrame
    tetramer_df = pd.DataFrame(tetramer_results)

    if not tetramer_df.empty:
        result_df = rescues_df.merge(
            tetramer_df,
            on="rescue_mutations",
            how="left"
        )
    else:
        result_df = rescues_df.copy()
        result_df["tetramer_impact_score"] = 0.0
        result_df["tetramer_risk"] = "unknown"

    # Filter high-risk if requested
    if filter_high_risk:
        high_risk_mask = result_df["tetramer_risk"].isin(["critical", "high"])
        n_filtered = high_risk_mask.sum()
        if n_filtered > 0:
            logger.warning(
                f"Filtered {n_filtered} rescues with critical/high tetramer risk"
            )
            result_df = result_df[~high_risk_mask].copy()

    return result_df


def summarize_tetramer_analysis(
    analysis_df: pd.DataFrame,
    output_path: Path | None = None,
) -> dict[str, Any]:
    """
    Generate summary statistics for tetramer interface analysis.

    Args:
        analysis_df: DataFrame from analyze_rescue_tetramer_impact()
        output_path: Optional path to save JSON summary

    Returns:
        Dict with summary statistics
    """
    summary = {
        "n_rescues_analyzed": len(analysis_df),
        "n_at_interface": int(analysis_df.get("at_interface_any", pd.Series([False])).sum()),
        "risk_distribution": {},
        "mean_impact_score": float(analysis_df.get("tetramer_impact_score", pd.Series([0])).mean()),
        "max_impact_score": float(analysis_df.get("tetramer_impact_score", pd.Series([0])).max()),
    }

    # Risk distribution
    if "tetramer_risk" in analysis_df.columns:
        risk_counts = analysis_df["tetramer_risk"].value_counts().to_dict()
        summary["risk_distribution"] = {str(k): int(v) for k, v in risk_counts.items()}

    # Top safe rescues
    safe_rescues = analysis_df[
        analysis_df["tetramer_risk"].isin(["safe", "low", "medium"])
    ].head(10)

    summary["top_safe_rescues"] = safe_rescues[[
        "rescue_mutations",
        "tetramer_impact_score",
        "tetramer_risk",
    ]].to_dict("records") if not safe_rescues.empty else []

    # Save if requested
    if output_path:
        with open(output_path, "w") as f:
            json.dump(summary, f, indent=2)
        logger.info(f"Saved tetramer summary to {output_path}")

    return summary
