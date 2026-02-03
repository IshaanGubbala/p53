"""
Multi-structure consensus scoring for robust ΔΔG predictions.

STRUCTURE WEIGHTING GUIDELINES:
-------------------------------
Experimental structures should be weighted HIGHER than computational predictions
because they represent validated conformations. Recommended weights:

Structure Type          | Resolution  | Recommended Weight
------------------------|-------------|-------------------
X-ray crystallography   | < 2.0 Å     | 3.0
X-ray crystallography   | 2.0-2.5 Å   | 2.5
X-ray crystallography   | 2.5-3.0 Å   | 2.0
Cryo-EM                 | < 3.0 Å     | 2.0
Cryo-EM                 | 3.0-4.0 Å   | 1.5
NMR ensemble (first)    | n/a         | 1.5
AlphaFold2 (high conf)  | pLDDT > 90  | 1.0
AlphaFold2 (med conf)   | pLDDT 70-90 | 0.7
AlphaFold2 (low conf)   | pLDDT < 70  | 0.3
ESMFold                 | n/a         | 0.5

DISCORDANT PREDICTIONS:
-----------------------
When structures disagree significantly (>2 kcal/mol), this indicates:
- Conformational flexibility at that site
- Possible artifacts in one structure
- Need for ensemble-based scoring

The code flags discordant predictions and can optionally exclude outliers.
"""
from __future__ import annotations

import logging
import warnings
from pathlib import Path
from typing import Optional, Callable
import numpy as np

from .evoef2_runner import score_mutation_set


logger = logging.getLogger(__name__)


# Default weights by structure type
DEFAULT_STRUCTURE_WEIGHTS = {
    # Experimental (high confidence)
    "xray_high": 3.0,      # X-ray < 2.0 Å
    "xray_medium": 2.5,    # X-ray 2.0-2.5 Å
    "xray_low": 2.0,       # X-ray 2.5-3.0 Å
    "cryoem_high": 2.0,    # Cryo-EM < 3.0 Å
    "cryoem_medium": 1.5,  # Cryo-EM 3.0-4.0 Å
    "nmr": 1.5,            # NMR ensemble

    # Computational (lower confidence)
    "alphafold_high": 1.0,   # pLDDT > 90
    "alphafold_medium": 0.7, # pLDDT 70-90
    "alphafold_low": 0.3,    # pLDDT < 70
    "esmfold": 0.5,
    "modeller": 0.5,

    # Default for unknown
    "experimental": 2.0,
    "predicted": 0.7,
    "default": 1.0,
}


def suggest_weight(
    structure_type: str,
    resolution: Optional[float] = None,
    plddt: Optional[float] = None,
) -> float:
    """
    Suggest appropriate weight for a structure based on type and quality.

    Args:
        structure_type: "xray", "cryoem", "nmr", "alphafold", "esmfold", etc.
        resolution: Resolution in Angstroms (for experimental structures)
        plddt: Mean pLDDT score (for AlphaFold/ESMFold)

    Returns:
        Suggested weight (higher = more trustworthy)
    """
    struct_type_lower = structure_type.lower()

    if "xray" in struct_type_lower or "crystal" in struct_type_lower:
        if resolution is not None:
            if resolution < 2.0:
                return DEFAULT_STRUCTURE_WEIGHTS["xray_high"]
            elif resolution < 2.5:
                return DEFAULT_STRUCTURE_WEIGHTS["xray_medium"]
            else:
                return DEFAULT_STRUCTURE_WEIGHTS["xray_low"]
        return DEFAULT_STRUCTURE_WEIGHTS["experimental"]

    elif "cryoem" in struct_type_lower or "cryo" in struct_type_lower:
        if resolution is not None:
            if resolution < 3.0:
                return DEFAULT_STRUCTURE_WEIGHTS["cryoem_high"]
            else:
                return DEFAULT_STRUCTURE_WEIGHTS["cryoem_medium"]
        return DEFAULT_STRUCTURE_WEIGHTS["cryoem_medium"]

    elif "nmr" in struct_type_lower:
        return DEFAULT_STRUCTURE_WEIGHTS["nmr"]

    elif "alphafold" in struct_type_lower or "af2" in struct_type_lower:
        if plddt is not None:
            if plddt > 90:
                return DEFAULT_STRUCTURE_WEIGHTS["alphafold_high"]
            elif plddt > 70:
                return DEFAULT_STRUCTURE_WEIGHTS["alphafold_medium"]
            else:
                return DEFAULT_STRUCTURE_WEIGHTS["alphafold_low"]
        return DEFAULT_STRUCTURE_WEIGHTS["alphafold_medium"]

    elif "esmfold" in struct_type_lower or "esm" in struct_type_lower:
        return DEFAULT_STRUCTURE_WEIGHTS["esmfold"]

    elif "modeller" in struct_type_lower or "homology" in struct_type_lower:
        return DEFAULT_STRUCTURE_WEIGHTS["modeller"]

    # Try to infer from name
    if any(x in struct_type_lower for x in ["pdb", "exp", "2ocj", "1tsr", "3kmd"]):
        return DEFAULT_STRUCTURE_WEIGHTS["experimental"]
    elif any(x in struct_type_lower for x in ["predicted", "model", "computed"]):
        return DEFAULT_STRUCTURE_WEIGHTS["predicted"]

    return DEFAULT_STRUCTURE_WEIGHTS["default"]


def check_discordance(
    scores: list[float],
    structure_ids: list[str],
    threshold: float = 2.0,
) -> dict:
    """
    Check for discordant predictions between structures.

    Discordant predictions (>2 kcal/mol difference) indicate:
    - Conformational flexibility at that site
    - Possible artifacts in one structure
    - Need for more careful analysis

    Args:
        scores: ΔΔG scores from each structure
        structure_ids: Structure identifiers
        threshold: Discordance threshold in kcal/mol

    Returns:
        Dictionary with discordance analysis
    """
    if len(scores) < 2:
        return {"discordant": False}

    scores_arr = np.array(scores)
    score_range = float(np.max(scores_arr) - np.min(scores_arr))

    result = {
        "discordant": score_range > threshold,
        "range": score_range,
        "threshold": threshold,
        "min_score": float(np.min(scores_arr)),
        "max_score": float(np.max(scores_arr)),
        "min_structure": structure_ids[int(np.argmin(scores_arr))],
        "max_structure": structure_ids[int(np.argmax(scores_arr))],
    }

    if result["discordant"]:
        result["warning"] = (
            f"Discordant predictions: {result['min_structure']} gives "
            f"ΔΔG={result['min_score']:.2f} but {result['max_structure']} gives "
            f"ΔΔG={result['max_score']:.2f} (difference={score_range:.2f} kcal/mol)"
        )
        warnings.warn(result["warning"], UserWarning)

    return result


def median_consensus(scores: list[float]) -> float:
    """Compute median of scores (robust to outliers).

    Args:
        scores: List of ΔΔG scores from different structures

    Returns:
        Median ΔΔG
    """
    return float(np.median(scores))


def mean_consensus(scores: list[float]) -> float:
    """Compute mean of scores.

    Args:
        scores: List of ΔΔG scores from different structures

    Returns:
        Mean ΔΔG
    """
    return float(np.mean(scores))


def weighted_consensus(scores: list[float], weights: list[float]) -> float:
    """Compute weighted average of scores.

    Args:
        scores: List of ΔΔG scores
        weights: List of weights (e.g., based on resolution/quality)

    Returns:
        Weighted mean ΔΔG
    """
    if len(scores) != len(weights):
        raise ValueError("scores and weights must have same length")

    total_weight = sum(weights)
    if total_weight == 0:
        return float(np.mean(scores))

    weighted_sum = sum(s * w for s, w in zip(scores, weights))
    return weighted_sum / total_weight


def score_mutation_set_multi(
    mutations: list[str],
    structures: list[dict],
    evoef2_cfg: dict,
    cache_dir: Path,
    work_root: Path,
    consensus_method: str = "weighted",  # Changed default to weighted
    require_all: bool = False,
    discordance_threshold: float = 2.0,
    auto_weight: bool = True,
) -> dict[str, any]:
    """Score mutations on multiple structures and compute consensus.

    IMPORTANT: Experimental structures should be weighted higher than predictions!
    Use auto_weight=True (default) to automatically assign weights based on
    structure type, or set "weight" in each structure config manually.

    Args:
        mutations: List of mutation tokens (e.g., ["A123V", "S95T"])
        structures: List of structure configs, each with keys:
            - id: Structure identifier (e.g., "alphafold", "2ocj_core")
            - pdb: Path to PDB file
            - chain_id: Chain identifier
            - weight: Weight for weighted consensus (optional, auto-assigned if not set)
            - type: Structure type for auto-weighting ("xray", "alphafold", etc.)
            - resolution: Resolution in Å (for experimental structures)
            - plddt: Mean pLDDT (for AlphaFold/ESMFold)
        evoef2_cfg: EvoEF2 configuration
        cache_dir: Cache directory
        work_root: Working directory root
        consensus_method: Consensus method ("median", "mean", "weighted")
            - "weighted" (RECOMMENDED): Uses structure-quality weights
            - "median": Robust to outliers, ignores weights
            - "mean": Simple average, ignores weights
        require_all: If True, fail if any structure fails to score
        discordance_threshold: Flag predictions differing by more than this (kcal/mol)
        auto_weight: Automatically assign weights based on structure type

    Returns:
        Dictionary with:
            - ddg_consensus: Consensus ΔΔG
            - ddg_<structure_id>: Per-structure ΔΔG scores
            - structures_scored: Number of successful scores
            - structures_failed: List of failed structure IDs
            - discordance: Discordance analysis (if multiple structures)
            - weights_used: Weights used for each structure
    """
    logger.info(f"Scoring {len(mutations)} mutations on {len(structures)} structures")

    results = {}
    scores = []
    weights = []
    structure_ids = []
    failed_structures = []
    weights_used = {}

    for struct_cfg in structures:
        struct_id = struct_cfg["id"]
        pdb_path = Path(struct_cfg["pdb"])
        chain_id = struct_cfg.get("chain_id", "A")

        # Determine weight
        if "weight" in struct_cfg:
            weight = struct_cfg["weight"]
        elif auto_weight:
            # Auto-assign weight based on structure type
            struct_type = struct_cfg.get("type", struct_id)
            resolution = struct_cfg.get("resolution")
            plddt = struct_cfg.get("plddt")
            weight = suggest_weight(struct_type, resolution, plddt)
            logger.info(f"    Auto-assigned weight {weight:.2f} for {struct_id} (type={struct_type})")
        else:
            weight = 1.0

        weights_used[struct_id] = weight

        logger.info(f"  Scoring on {struct_id} ({pdb_path.name}, chain {chain_id}, weight={weight:.2f})")

        try:
            # Score on this structure
            ddg = score_mutation_set(
                mutations=mutations,
                pdb_path=pdb_path,
                cache_dir=cache_dir,
                evoef2_cfg=evoef2_cfg,
                work_root=work_root,
            )

            results[f"ddg_{struct_id}"] = ddg
            scores.append(ddg)
            weights.append(weight)
            structure_ids.append(struct_id)

            logger.info(f"    ΔΔG = {ddg:.3f}")

        except Exception as e:
            logger.warning(f"    Failed to score on {struct_id}: {e}")
            results[f"ddg_{struct_id}"] = None
            failed_structures.append(struct_id)

            if require_all:
                raise RuntimeError(f"Scoring failed on {struct_id} (require_all=True)")

    # Check if we have any successful scores
    if not scores:
        raise RuntimeError(f"All structures failed to score (tried {len(structures)})")

    logger.info(f"  Scored on {len(scores)}/{len(structures)} structures")

    # Compute consensus
    if consensus_method == "median":
        consensus = median_consensus(scores)
    elif consensus_method == "mean":
        consensus = mean_consensus(scores)
    elif consensus_method == "weighted":
        consensus = weighted_consensus(scores, weights)
    else:
        logger.warning(f"Unknown consensus method '{consensus_method}', using median")
        consensus = median_consensus(scores)

    logger.info(f"  Consensus ΔΔG ({consensus_method}): {consensus:.3f}")

    # Add consensus and metadata
    results["ddg_consensus"] = consensus
    results["structures_scored"] = len(scores)
    results["structures_failed"] = failed_structures
    results["consensus_method"] = consensus_method
    results["weights_used"] = weights_used

    # Add score statistics
    if len(scores) > 1:
        results["ddg_std"] = float(np.std(scores))
        results["ddg_min"] = float(np.min(scores))
        results["ddg_max"] = float(np.max(scores))
        results["ddg_range"] = float(np.max(scores) - np.min(scores))

        # Check for discordant predictions
        discordance = check_discordance(scores, structure_ids, discordance_threshold)
        results["discordance"] = discordance

        if discordance["discordant"]:
            logger.warning(f"  ⚠️ {discordance['warning']}")

    return results


def validate_structures(
    structures: list[dict],
    check_compatibility: bool = True,
) -> dict[str, any]:
    """Validate structure configurations before scoring.

    Args:
        structures: List of structure configs
        check_compatibility: If True, check structures are compatible

    Returns:
        Dictionary with validation results
    """
    from .structure_validator import get_structure_info, compare_structures, validate_residue_numbering

    validation = {
        "valid": True,
        "structures": {},
        "warnings": [],
    }

    # Check each structure exists and is valid
    for struct_cfg in structures:
        struct_id = struct_cfg["id"]
        pdb_path = Path(struct_cfg["pdb"])
        chain_id = struct_cfg.get("chain_id", "A")

        if not pdb_path.exists():
            validation["valid"] = False
            validation["warnings"].append(f"{struct_id}: PDB file not found at {pdb_path}")
            continue

        info = get_structure_info(pdb_path, chain_id)
        validation["structures"][struct_id] = info

        if "error" in info:
            validation["valid"] = False
            validation["warnings"].append(f"{struct_id}: {info['error']}")

    # Check compatibility between structures
    if check_compatibility and len(structures) > 1:
        ref_struct = structures[0]
        ref_id = ref_struct["id"]
        ref_pdb = Path(ref_struct["pdb"])
        ref_chain = ref_struct.get("chain_id", "A")

        for struct_cfg in structures[1:]:
            struct_id = struct_cfg["id"]
            pdb_path = Path(struct_cfg["pdb"])
            chain_id = struct_cfg.get("chain_id", "A")

            comparison = compare_structures(ref_pdb, pdb_path, ref_chain, chain_id)

            if not comparison.get("compatible", False):
                validation["warnings"].append(
                    f"{ref_id} vs {struct_id}: Low sequence identity "
                    f"({comparison.get('sequence_identity', 0):.1%})"
                )

            validation[f"comparison_{ref_id}_vs_{struct_id}"] = comparison

    return validation
