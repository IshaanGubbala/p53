"""
Druggability analysis for rescue mutations.

Alternative to fpocket that uses BioPython + SASA data to identify
potential druggable pockets near rescue mutations.
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


def calculate_pocket_score(
    pdb_path: Path,
    rescue_positions: list[int],
    sasa_data: dict[str, Any] | None = None,
    distance_threshold: float = 12.0,
) -> dict[str, Any]:
    """
    Calculate druggability score for regions near rescue mutations.

    A "druggable pocket" is characterized by:
    1. Surface accessibility (SASA > 20 Å²)
    2. Hydrophobic character (aromatic/aliphatic residues)
    3. Concavity (surrounded by other residues)
    4. Size (enough space for small molecule binding)

    Args:
        pdb_path: Path to PDB structure
        rescue_positions: List of residue numbers for rescue mutations
        sasa_data: Pre-computed SASA data (optional)
        distance_threshold: Radius around rescue site to analyze (Å)

    Returns:
        Dictionary with pocket metrics
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", str(pdb_path))
    model = structure[0]

    # Get all atoms for neighbor search
    atoms = [atom for atom in model.get_atoms()]
    ns = NeighborSearch(atoms)

    results = {}

    for rescue_pos in rescue_positions:
        try:
            # Find rescue residue
            rescue_residue = None
            for chain in model:
                for residue in chain:
                    if residue.id[1] == rescue_pos:
                        rescue_residue = residue
                        break
                if rescue_residue:
                    break

            if not rescue_residue:
                logger.warning(f"Rescue position {rescue_pos} not found in structure")
                continue

            # Get CA atom for center point
            if "CA" not in rescue_residue:
                logger.warning(f"No CA atom for position {rescue_pos}")
                continue

            center_atom = rescue_residue["CA"]

            # Find nearby residues
            nearby_residues = set()
            for atom in ns.search(center_atom.coord, distance_threshold):
                residue = atom.get_parent()
                if residue.id[0] == " ":  # Standard residue
                    nearby_residues.add(residue)

            # Analyze pocket characteristics
            pocket_residues = []
            hydrophobic_count = 0
            aromatic_count = 0
            charged_count = 0
            polar_count = 0

            hydrophobic = {"ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TRP", "PRO"}
            aromatic = {"PHE", "TRP", "TYR", "HIS"}
            charged = {"ARG", "LYS", "ASP", "GLU", "HIS"}
            polar = {"SER", "THR", "ASN", "GLN", "CYS", "TYR"}

            for residue in nearby_residues:
                resname = residue.get_resname()
                res_id = residue.id[1]

                pocket_residues.append({
                    "position": res_id,
                    "residue": resname,
                    "distance": (residue["CA"].coord - center_atom.coord).sum() if "CA" in residue else None
                })

                if resname in hydrophobic:
                    hydrophobic_count += 1
                if resname in aromatic:
                    aromatic_count += 1
                if resname in charged:
                    charged_count += 1
                if resname in polar:
                    polar_count += 1

            # Calculate druggability score
            # Higher score = more druggable
            n_residues = len(pocket_residues)

            if n_residues == 0:
                druggability_score = 0.0
            else:
                # Score components (normalized 0-1)
                size_score = min(n_residues / 20.0, 1.0)  # Optimal ~15-20 residues
                hydrophobic_ratio = hydrophobic_count / n_residues
                aromatic_ratio = aromatic_count / n_residues

                # Druggability formula (empirical weights)
                druggability_score = (
                    0.4 * size_score +
                    0.3 * hydrophobic_ratio +
                    0.2 * aromatic_ratio +
                    0.1 * (1.0 - charged_count / max(n_residues, 1))  # Prefer uncharged
                )

            # Get SASA if available
            sasa_value = None
            if sasa_data and str(rescue_pos) in sasa_data:
                sasa_value = sasa_data[str(rescue_pos)]

            # Classify druggability
            if druggability_score >= 0.6:
                classification = "High"
                flag = "🎯 Druggable Pocket"
            elif druggability_score >= 0.4:
                classification = "Medium"
                flag = "⚠️ Moderate Pocket"
            else:
                classification = "Low"
                flag = "❌ Poor Druggability"

            results[rescue_pos] = {
                "position": rescue_pos,
                "druggability_score": round(druggability_score, 3),
                "classification": classification,
                "flag": flag,
                "n_pocket_residues": n_residues,
                "hydrophobic_count": hydrophobic_count,
                "aromatic_count": aromatic_count,
                "charged_count": charged_count,
                "polar_count": polar_count,
                "hydrophobic_ratio": round(hydrophobic_count / max(n_residues, 1), 3),
                "aromatic_ratio": round(aromatic_count / max(n_residues, 1), 3),
                "rescue_sasa": round(sasa_value, 2) if sasa_value else None,
                "pocket_residues": pocket_residues[:10],  # Top 10 for brevity
            }

        except Exception as e:
            logger.error(f"Error analyzing pocket at position {rescue_pos}: {e}")
            continue

    return results


def check_functional_sites(
    positions: list[int],
    functional_sites: dict[str, list[int]],
) -> dict[str, Any]:
    """
    Check if rescue positions overlap with functionally critical sites.

    Args:
        positions: List of residue positions to check
        functional_sites: Dict of functional site categories with residue lists

    Returns:
        Dictionary with functional site warnings and flags
    """
    warnings = []
    violated_categories = []
    all_functional = set()

    # Flatten all functional sites
    for category, residues in functional_sites.items():
        all_functional.update(residues)

    # Check each position
    for pos in positions:
        if pos in all_functional:
            # Find which categories this position belongs to
            for category, residues in functional_sites.items():
                if pos in residues:
                    warnings.append(f"{pos} ({category})")
                    violated_categories.append(category)

    return {
        "has_functional_site_conflict": len(warnings) > 0,
        "functional_site_warnings": "; ".join(warnings) if warnings else None,
        "violated_categories": list(set(violated_categories)),
        "n_functional_conflicts": len(warnings),
    }


def analyze_rescue_druggability(
    rescues_df: pd.DataFrame,
    target: str,
    base_pdb: Path,
    sasa_path: Path | None = None,
    top_n: int = 10,
    evoef2_cfg: dict | None = None,
    work_dir: Path | None = None,
    functional_sites: dict[str, list[int]] | None = None,
    filter_functional: bool = True,
) -> pd.DataFrame:
    """
    Analyze druggability of top rescue mutations.

    IMPORTANT: Builds rescued structures (cancer + rescue) using EvoEF2 BuildMutant
    to analyze pockets in the RESCUED state, not just the cancer mutant.

    NEW: Filters out rescues that mutate functionally critical residues (DNA contacts,
    phosphorylation sites, zinc coordination, etc.) to avoid "Pyrrhic victories" where
    we create a druggable pocket by destroying protein function.

    Args:
        rescues_df: DataFrame with rescue candidates
        target: Target mutation name (e.g., "R175H")
        base_pdb: Path to wild-type structure (used as template for building)
        sasa_path: Path to SASA JSON file
        top_n: Number of top rescues to analyze
        evoef2_cfg: EvoEF2 configuration for BuildMutant
        work_dir: Working directory for structure building
        functional_sites: Dict of functional site categories with residue lists
        filter_functional: If True, exclude rescues that hit functional sites

    Returns:
        DataFrame with added druggability columns and functional site warnings
    """
    from src.scoring.evoef2_runner import build_mutant_model

    # Load SASA data
    sasa_data = None
    if sasa_path and sasa_path.exists():
        with open(sasa_path) as f:
            sasa_data = json.load(f)

    # Setup for structure building
    if evoef2_cfg is None:
        evoef2_cfg = {}
    if work_dir is None:
        work_dir = Path("Data/processed/cache/druggability_structures")
    work_dir.mkdir(parents=True, exist_ok=True)

    # Get top rescues
    top_rescues = rescues_df.head(top_n).copy()

    # Extract rescue positions
    druggability_results = []

    for idx, row in top_rescues.iterrows():
        rescue_mutations = row.get("rescue_mutations", "")
        full_mutations = row.get("full_mutations", "")

        if pd.isna(rescue_mutations) or pd.isna(full_mutations):
            continue

        # Parse rescue positions (e.g., "M133L,A189S" -> [133, 189])
        positions = []
        for mut in str(rescue_mutations).split(","):
            mut = mut.strip()
            if len(mut) >= 2:
                try:
                    # Extract position number (e.g., M133L -> 133)
                    pos = int(''.join(filter(str.isdigit, mut)))
                    positions.append(pos)
                except ValueError:
                    continue

        if not positions:
            continue

        # Check for functional site conflicts
        functional_check = {"has_functional_site_conflict": False}
        if functional_sites:
            functional_check = check_functional_sites(positions, functional_sites)

            # Skip this rescue if filtering is enabled and it hits functional sites
            if filter_functional and functional_check["has_functional_site_conflict"]:
                logger.warning(
                    f"Skipping {rescue_mutations}: conflicts with functional sites "
                    f"({functional_check['functional_site_warnings']})"
                )
                continue

        # Build RESCUED structure (cancer + rescue mutations)
        try:
            full_muts = [m.strip() for m in str(full_mutations).split(",") if m.strip()]

            # Build rescued structure using EvoEF2
            rescued_pdb = build_mutant_model(
                full_muts,
                base_pdb,
                evoef2_cfg,
                work_dir,
                recompute=False,  # Use cache
            )

            logger.info(f"Built rescued structure for {rescue_mutations}: {rescued_pdb}")

        except Exception as e:
            logger.warning(f"Failed to build rescued structure for {rescue_mutations}: {e}")
            # Fallback to base structure (will give incorrect results but won't crash)
            rescued_pdb = base_pdb

        # Analyze pocket in RESCUED structure
        try:
            pocket_analysis = calculate_pocket_score(
                rescued_pdb,  # Use rescued structure, not base!
                positions,
                sasa_data,
                distance_threshold=12.0
            )

            # Aggregate scores for multi-mutation rescues
            if pocket_analysis:
                avg_score = np.mean([p["druggability_score"] for p in pocket_analysis.values()])
                max_score = max([p["druggability_score"] for p in pocket_analysis.values()])
                flags = [p["flag"] for p in pocket_analysis.values()]

                druggability_results.append({
                    "rescue_mutations": rescue_mutations,
                    "druggability_score": round(avg_score, 3),
                    "max_druggability": round(max_score, 3),
                    "pocket_flag": flags[0] if len(flags) == 1 else "🔀 Multi-site",
                    "n_positions_analyzed": len(pocket_analysis),
                    "pocket_details": pocket_analysis,
                    "functional_site_conflict": functional_check.get("has_functional_site_conflict", False),
                    "functional_warnings": functional_check.get("functional_site_warnings", None),
                    "violated_categories": functional_check.get("violated_categories", []),
                })

        except Exception as e:
            logger.warning(f"Failed to analyze druggability for {rescue_mutations}: {e}")
            druggability_results.append({
                "rescue_mutations": rescue_mutations,
                "druggability_score": 0.0,
                "max_druggability": 0.0,
                "pocket_flag": "❓ Analysis Failed",
                "n_positions_analyzed": 0,
                "pocket_details": {},
            })

    # Merge with original DataFrame
    druggability_df = pd.DataFrame(druggability_results)

    if not druggability_df.empty:
        result_df = top_rescues.merge(
            druggability_df,
            on="rescue_mutations",
            how="left"
        )
    else:
        result_df = top_rescues.copy()
        result_df["druggability_score"] = 0.0
        result_df["pocket_flag"] = "❓ Not Analyzed"

    return result_df


def summarize_druggability(
    druggability_df: pd.DataFrame,
    output_path: Path | None = None,
) -> dict[str, Any]:
    """
    Generate summary statistics for druggability analysis.

    Args:
        druggability_df: DataFrame with druggability scores
        output_path: Optional path to save JSON summary

    Returns:
        Dictionary with summary statistics
    """
    summary = {
        "n_rescues_analyzed": int(len(druggability_df)),
        "high_druggability": int((druggability_df["druggability_score"] >= 0.6).sum()),
        "medium_druggability": int(
            ((druggability_df["druggability_score"] >= 0.4) &
             (druggability_df["druggability_score"] < 0.6)).sum()
        ),
        "low_druggability": int((druggability_df["druggability_score"] < 0.4).sum()),
        "mean_score": round(float(druggability_df["druggability_score"].mean()), 3),
        "max_score": round(float(druggability_df["druggability_score"].max()), 3),
        "top_druggable_rescues": druggability_df.nlargest(3, "druggability_score")[
            ["rescue_mutations", "druggability_score", "pocket_flag"]
        ].to_dict("records"),
    }

    if output_path:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w") as f:
            json.dump(summary, f, indent=2)
        logger.info(f"Druggability summary saved to {output_path}")

    return summary
