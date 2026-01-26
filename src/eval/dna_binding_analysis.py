"""
DNA binding affinity analysis for rescue mutations.

This module identifies DNA-contacting residues and estimates the impact
of mutations on protein-DNA binding affinity.

CRITICAL: A rescue that stabilizes the protein but destroys DNA binding
is functionally worse than no rescue at all.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, NeighborSearch, Selection

from src.core.logging import get_logger

logger = get_logger(__name__)


def identify_dna_contacts(
    pdb_path: Path,
    protein_chains: list[str] = ["A"],
    dna_chains: list[str] = ["E", "F"],
    distance_threshold: float = 5.0,
) -> dict[int, dict[str, Any]]:
    """
    Identify residues in direct contact with DNA.

    Args:
        pdb_path: Path to PDB structure with protein-DNA complex
        protein_chains: List of protein chain IDs to analyze
        dna_chains: List of DNA chain IDs
        distance_threshold: Maximum distance (Å) for contact (default: 5.0)

    Returns:
        Dictionary mapping residue position to contact information
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", str(pdb_path))
    model = structure[0]

    # Get protein and DNA atoms
    protein_atoms = []
    dna_atoms = []

    for chain in model:
        chain_id = chain.get_id()
        if chain_id in protein_chains:
            for residue in chain:
                if residue.get_id()[0] == " ":  # Standard residue
                    protein_atoms.extend(residue.get_atoms())
        elif chain_id in dna_chains:
            for residue in chain:
                dna_atoms.extend(residue.get_atoms())

    logger.info(f"Found {len(protein_atoms)} protein atoms, {len(dna_atoms)} DNA atoms")

    # Build neighbor search for DNA
    dna_search = NeighborSearch(dna_atoms)

    # Find contacts
    contacts = {}

    for atom in protein_atoms:
        residue = atom.get_parent()
        res_id = residue.get_id()[1]  # Residue number
        res_name = residue.get_resname()

        # Find DNA atoms within threshold
        nearby_dna = dna_search.search(atom.coord, distance_threshold)

        if nearby_dna:
            if res_id not in contacts:
                contacts[res_id] = {
                    "position": int(res_id),
                    "residue_name": str(res_name),
                    "n_contacts": 0,
                    "min_distance": float('inf'),
                    "contacting_atoms": [],
                    "contact_types": set(),
                }

            # Calculate distances
            for dna_atom in nearby_dna:
                distance = float(np.linalg.norm(atom.coord - dna_atom.coord))
                if distance <= distance_threshold:
                    contacts[res_id]["n_contacts"] += 1
                    contacts[res_id]["min_distance"] = float(min(
                        contacts[res_id]["min_distance"], distance
                    ))
                    contacts[res_id]["contacting_atoms"].append({
                        "protein_atom": atom.get_name(),
                        "dna_atom": dna_atom.get_name(),
                        "distance": distance,
                    })

                    # Classify contact type
                    if atom.get_name() in ["N", "O"] or dna_atom.get_name() in ["N", "O"]:
                        contacts[res_id]["contact_types"].add("hydrogen_bond")
                    elif "O" in dna_atom.get_parent().get_resname():  # Phosphate
                        contacts[res_id]["contact_types"].add("electrostatic")
                    else:
                        contacts[res_id]["contact_types"].add("van_der_waals")

    # Convert sets to lists and ensure all types are JSON serializable
    for pos in contacts:
        contacts[pos]["contact_types"] = list(contacts[pos]["contact_types"])
        # Ensure min_distance is a regular float, not inf
        if contacts[pos]["min_distance"] == float('inf'):
            contacts[pos]["min_distance"] = None
        else:
            contacts[pos]["min_distance"] = float(contacts[pos]["min_distance"])

    logger.info(f"Identified {len(contacts)} DNA-contacting residues")

    return contacts


def classify_dna_contact_strength(
    contacts: dict[int, dict[str, Any]]
) -> dict[int, str]:
    """
    Classify DNA contacts as strong/medium/weak based on contact characteristics.

    Args:
        contacts: Dict from identify_dna_contacts()

    Returns:
        Dict mapping position to classification (strong/medium/weak)
    """
    classifications = {}

    for pos, contact_info in contacts.items():
        n_contacts = contact_info["n_contacts"]
        min_dist = contact_info["min_distance"]
        contact_types = contact_info["contact_types"]

        # Strong contact: many contacts, close distance, H-bonds or electrostatic
        if n_contacts >= 5 and min_dist < 3.5:
            classifications[pos] = "strong"
        elif n_contacts >= 3 and min_dist < 4.0 and (
            "hydrogen_bond" in contact_types or "electrostatic" in contact_types
        ):
            classifications[pos] = "strong"
        # Medium contact: moderate contacts
        elif n_contacts >= 2 and min_dist < 4.5:
            classifications[pos] = "medium"
        # Weak contact: few contacts or large distance
        else:
            classifications[pos] = "weak"

    return classifications


def estimate_binding_impact(
    mutation: str,
    dna_contacts: dict[int, dict[str, Any]],
    contact_classifications: dict[int, str],
) -> dict[str, Any]:
    """
    Estimate the impact of a mutation on DNA binding.

    This is a heuristic approximation. For rigorous analysis,
    use FoldX, Rosetta, or HADDOCK for actual ΔΔG calculations.

    Args:
        mutation: Mutation string (e.g., "R196Q")
        dna_contacts: DNA contact information
        contact_classifications: Contact strength classifications

    Returns:
        Dict with impact assessment
    """
    # Parse mutation
    if len(mutation) < 3:
        return {"error": "Invalid mutation format"}

    wt_aa = mutation[0]
    position = int(''.join(filter(str.isdigit, mutation)))
    mut_aa = mutation[-1]

    # Check if position contacts DNA
    if position not in dna_contacts:
        return {
            "position": int(position),
            "contacts_dna": False,
            "impact": "none",
            "impact_score": 0.0,
            "risk_level": "safe",
            "explanation": f"Position {position} does not contact DNA",
        }

    # Get contact info
    contact_info = dna_contacts[position]
    contact_strength = contact_classifications.get(position, "unknown")

    # Amino acid properties
    charged_positive = ["R", "K", "H"]
    charged_negative = ["D", "E"]
    polar = ["S", "T", "N", "Q"]
    hydrophobic = ["A", "V", "L", "I", "M", "F", "W", "Y"]

    # Assess mutation impact
    impact_score = 0
    explanations = []

    # 1. Loss of positive charge (critical for DNA phosphate backbone)
    if wt_aa in charged_positive and mut_aa not in charged_positive:
        impact_score += 3
        explanations.append(
            f"Loss of positive charge ({wt_aa}→{mut_aa}) weakens electrostatic DNA interaction"
        )

    # 2. Loss of H-bond donor/acceptor
    if "hydrogen_bond" in contact_info["contact_types"]:
        if (wt_aa in polar or wt_aa in charged_positive + charged_negative) and \
           mut_aa not in (polar + charged_positive + charged_negative):
            impact_score += 2
            explanations.append(f"Loss of H-bond capability ({wt_aa}→{mut_aa})")

    # 3. Introduction of steric clash or unfavorable charge
    if mut_aa in charged_negative and "electrostatic" in contact_info["contact_types"]:
        impact_score += 2
        explanations.append(f"Introduction of negative charge near DNA phosphate (repulsion)")

    # 4. Contact strength multiplier
    if contact_strength == "strong":
        impact_score *= 1.5
    elif contact_strength == "medium":
        impact_score *= 1.2

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
        "contacts_dna": True,
        "contact_strength": str(contact_strength),
        "n_contacts": int(contact_info["n_contacts"]),
        "min_distance": float(contact_info["min_distance"]) if contact_info["min_distance"] is not None else None,
        "contact_types": list(contact_info["contact_types"]),
        "impact_score": float(round(impact_score, 2)),
        "impact": str(impact),
        "risk_level": str(risk_level),
        "explanations": str("; ".join(explanations) if explanations else "Minimal impact expected"),
    }


def analyze_rescue_dna_binding(
    rescues_df: pd.DataFrame,
    dna_complex_pdb: Path,
    protein_chains: list[str] = ["A"],
    dna_chains: list[str] = ["E", "F"],
    filter_high_risk: bool = True,
) -> pd.DataFrame:
    """
    Analyze DNA binding impact for all rescue mutations.

    Args:
        rescues_df: DataFrame with rescue candidates
        dna_complex_pdb: Path to protein-DNA complex structure
        protein_chains: Protein chain IDs
        dna_chains: DNA chain IDs
        filter_high_risk: If True, flag rescues with high DNA binding risk

    Returns:
        DataFrame with DNA binding impact columns
    """
    # Identify DNA contacts
    logger.info("Identifying DNA-contacting residues...")
    dna_contacts = identify_dna_contacts(
        dna_complex_pdb,
        protein_chains=protein_chains,
        dna_chains=dna_chains,
    )

    # Classify contact strengths
    contact_classifications = classify_dna_contact_strength(dna_contacts)

    logger.info(f"Found {len(dna_contacts)} DNA-contacting residues:")
    for pos in sorted(dna_contacts.keys()):
        strength = contact_classifications.get(pos, "unknown")
        n_contacts = dna_contacts[pos]["n_contacts"]
        logger.info(f"  Position {pos}: {strength} contact ({n_contacts} interactions)")

    # Analyze each rescue
    binding_results = []

    for idx, row in rescues_df.iterrows():
        rescue_mutations = row.get("rescue_mutations", "")

        if pd.isna(rescue_mutations):
            continue

        # Parse mutations
        mutations = [m.strip() for m in str(rescue_mutations).split(",")]

        # Check each mutation for DNA binding impact
        max_impact_score = 0
        max_risk_level = "safe"
        impacts = []

        for mut in mutations:
            if len(mut) < 3:
                continue

            impact = estimate_binding_impact(
                mut,
                dna_contacts,
                contact_classifications,
            )

            impacts.append(impact)

            if impact.get("impact_score", 0) > max_impact_score:
                max_impact_score = impact["impact_score"]
                max_risk_level = impact["risk_level"]

        binding_results.append({
            "rescue_mutations": rescue_mutations,
            "dna_binding_impact_score": max_impact_score,
            "dna_binding_risk": max_risk_level,
            "dna_contacts_any": any(imp.get("contacts_dna", False) for imp in impacts),
            "dna_binding_details": impacts,
        })

    # Merge with original DataFrame
    binding_df = pd.DataFrame(binding_results)

    if not binding_df.empty:
        result_df = rescues_df.merge(
            binding_df,
            on="rescue_mutations",
            how="left"
        )
    else:
        result_df = rescues_df.copy()
        result_df["dna_binding_impact_score"] = 0.0
        result_df["dna_binding_risk"] = "unknown"

    # Filter high-risk if requested
    if filter_high_risk:
        high_risk_mask = result_df["dna_binding_risk"].isin(["critical", "high"])
        n_filtered = high_risk_mask.sum()
        if n_filtered > 0:
            logger.warning(
                f"Filtered {n_filtered} rescues with critical/high DNA binding risk"
            )
            result_df = result_df[~high_risk_mask].copy()

    return result_df


def summarize_dna_binding_analysis(
    analysis_df: pd.DataFrame,
    output_path: Path | None = None,
) -> dict[str, Any]:
    """
    Generate summary statistics for DNA binding analysis.

    Args:
        analysis_df: DataFrame from analyze_rescue_dna_binding()
        output_path: Optional path to save JSON summary

    Returns:
        Dict with summary statistics
    """
    summary = {
        "n_rescues_analyzed": len(analysis_df),
        "n_contact_dna": int(analysis_df.get("dna_contacts_any", pd.Series([False])).sum()),
        "risk_distribution": {},
        "mean_impact_score": float(analysis_df.get("dna_binding_impact_score", pd.Series([0])).mean()),
        "max_impact_score": float(analysis_df.get("dna_binding_impact_score", pd.Series([0])).max()),
    }

    # Risk distribution
    if "dna_binding_risk" in analysis_df.columns:
        risk_counts = analysis_df["dna_binding_risk"].value_counts().to_dict()
        summary["risk_distribution"] = {str(k): int(v) for k, v in risk_counts.items()}

    # Top safe rescues (low DNA binding risk)
    safe_rescues = analysis_df[
        analysis_df["dna_binding_risk"].isin(["safe", "low", "medium"])
    ].head(10)

    summary["top_safe_rescues"] = safe_rescues[[
        "rescue_mutations",
        "dna_binding_impact_score",
        "dna_binding_risk",
    ]].to_dict("records") if not safe_rescues.empty else []

    # Save if requested
    if output_path:
        with open(output_path, "w") as f:
            json.dump(summary, f, indent=2)
        logger.info(f"Saved DNA binding summary to {output_path}")

    return summary
