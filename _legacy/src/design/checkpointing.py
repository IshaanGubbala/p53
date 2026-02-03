"""Checkpointing utilities for beam search rescue design."""
from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

import pandas as pd

from src.core.hashing import dict_sha256


logger = logging.getLogger(__name__)


def compute_config_hash(config: dict[str, Any]) -> str:
    """Compute hash of configuration for checkpoint validation.

    Args:
        config: Configuration dictionary

    Returns:
        SHA256 hash of config
    """
    return dict_sha256(config)


def save_checkpoint(
    checkpoint_dir: Path,
    step: int,
    candidates: pd.DataFrame,
    beam_state: list[tuple[str, ...]],
    config_hash: str,
    metadata: dict[str, Any] | None = None,
) -> None:
    """Save checkpoint after beam search step.

    Args:
        checkpoint_dir: Directory to save checkpoints
        step: Current step number (1, 2, 3, ...)
        candidates: DataFrame with all scored candidates so far
        beam_state: Current beam (top candidates for expansion)
        config_hash: Hash of configuration for validation
        metadata: Optional metadata to save
    """
    checkpoint_dir.mkdir(parents=True, exist_ok=True)

    # Save candidates dataframe
    candidates_file = checkpoint_dir / f"step_{step}.parquet"
    candidates.to_parquet(candidates_file, index=False)

    # Save beam state and metadata
    meta_file = checkpoint_dir / "checkpoint_meta.json"
    checkpoint_meta = {
        "current_step": step,
        "config_hash": config_hash,
        "beam_size": len(beam_state),
        "n_candidates": len(candidates),
        "beam_state": [list(mset) for mset in beam_state],
    }

    if metadata:
        checkpoint_meta.update(metadata)

    meta_file.write_text(json.dumps(checkpoint_meta, indent=2), encoding="utf-8")

    logger.info(f"Saved checkpoint at step {step} ({len(candidates)} candidates, beam size {len(beam_state)})")


def load_checkpoint(
    checkpoint_dir: Path,
    expected_config_hash: str | None = None,
    verify_config: bool = True,
) -> dict[str, Any] | None:
    """Load checkpoint from directory.

    Args:
        checkpoint_dir: Directory containing checkpoints
        expected_config_hash: Expected config hash for validation
        verify_config: If True, verify config hash matches

    Returns:
        Dictionary with checkpoint data, or None if no valid checkpoint
    """
    if not checkpoint_dir.exists():
        return None

    meta_file = checkpoint_dir / "checkpoint_meta.json"
    if not meta_file.exists():
        return None

    # Load metadata
    meta = json.loads(meta_file.read_text(encoding="utf-8"))
    step = meta.get("current_step")

    if step is None:
        logger.warning("Checkpoint metadata missing 'current_step', ignoring checkpoint")
        return None

    # Verify config hash
    if verify_config and expected_config_hash:
        stored_hash = meta.get("config_hash")
        if stored_hash != expected_config_hash:
            logger.warning(
                f"Config hash mismatch: checkpoint={stored_hash[:8]}... vs current={expected_config_hash[:8]}..."
            )
            logger.warning("Configuration has changed, ignoring checkpoint")
            logger.warning("Use --recompute to force fresh run")
            return None

    # Load candidates
    candidates_file = checkpoint_dir / f"step_{step}.parquet"
    if not candidates_file.exists():
        logger.warning(f"Checkpoint step_{step}.parquet not found, ignoring checkpoint")
        return None

    candidates = pd.read_parquet(candidates_file)

    # Reconstruct beam state
    beam_state_raw = meta.get("beam_state", [])
    beam_state = [tuple(mset) for mset in beam_state_raw]

    logger.info(f"Loaded checkpoint from step {step} ({len(candidates)} candidates, beam size {len(beam_state)})")

    return {
        "step": step,
        "candidates": candidates,
        "beam_state": beam_state,
        "metadata": meta,
    }


def delete_checkpoint(checkpoint_dir: Path) -> None:
    """Delete all checkpoint files in directory.

    Args:
        checkpoint_dir: Directory containing checkpoints
    """
    if not checkpoint_dir.exists():
        return

    # Delete checkpoint files
    for file in checkpoint_dir.glob("step_*.parquet"):
        file.unlink()

    meta_file = checkpoint_dir / "checkpoint_meta.json"
    if meta_file.exists():
        meta_file.unlink()

    logger.info(f"Deleted checkpoints from {checkpoint_dir}")


def should_resume(
    checkpoint_dir: Path,
    config_hash: str,
    force_recompute: bool = False,
) -> bool:
    """Determine if checkpoint should be resumed.

    Args:
        checkpoint_dir: Directory containing checkpoints
        config_hash: Current config hash
        force_recompute: If True, always return False

    Returns:
        True if checkpoint exists and should be resumed
    """
    if force_recompute:
        return False

    checkpoint = load_checkpoint(checkpoint_dir, config_hash, verify_config=True)
    return checkpoint is not None


def get_checkpoint_step(checkpoint_dir: Path) -> int | None:
    """Get the current step from checkpoint metadata.

    Args:
        checkpoint_dir: Directory containing checkpoints

    Returns:
        Current step number, or None if no checkpoint
    """
    meta_file = checkpoint_dir / "checkpoint_meta.json"
    if not meta_file.exists():
        return None

    try:
        meta = json.loads(meta_file.read_text(encoding="utf-8"))
        return meta.get("current_step")
    except Exception:
        return None
