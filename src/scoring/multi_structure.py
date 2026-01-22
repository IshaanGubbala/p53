"""Multi-structure consensus scoring for robust ΔΔG predictions."""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional, Callable
import numpy as np

from .evoef2_runner import score_mutation_set


logger = logging.getLogger(__name__)


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
    consensus_method: str = "median",
    require_all: bool = False,
) -> dict[str, any]:
    """Score mutations on multiple structures and compute consensus.

    Args:
        mutations: List of mutation tokens (e.g., ["A123V", "S95T"])
        structures: List of structure configs, each with keys:
            - id: Structure identifier (e.g., "alphafold", "2ocj_core")
            - pdb: Path to PDB file
            - chain_id: Chain identifier
            - weight: Weight for weighted consensus (optional)
        evoef2_cfg: EvoEF2 configuration
        cache_dir: Cache directory
        work_root: Working directory root
        consensus_method: Consensus method ("median", "mean", "weighted")
        require_all: If True, fail if any structure fails to score

    Returns:
        Dictionary with:
            - ddg_consensus: Consensus ΔΔG
            - ddg_<structure_id>: Per-structure ΔΔG scores
            - structures_scored: Number of successful scores
            - structures_failed: List of failed structure IDs
    """
    logger.info(f"Scoring {len(mutations)} mutations on {len(structures)} structures")

    results = {}
    scores = []
    weights = []
    failed_structures = []

    for struct_cfg in structures:
        struct_id = struct_cfg["id"]
        pdb_path = Path(struct_cfg["pdb"])
        chain_id = struct_cfg.get("chain_id", "A")
        weight = struct_cfg.get("weight", 1.0)

        logger.info(f"  Scoring on {struct_id} ({pdb_path.name}, chain {chain_id})")

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

    # Add score statistics
    if len(scores) > 1:
        results["ddg_std"] = float(np.std(scores))
        results["ddg_min"] = float(np.min(scores))
        results["ddg_max"] = float(np.max(scores))
        results["ddg_range"] = float(np.max(scores) - np.min(scores))

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
