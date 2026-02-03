"""Cache utilities for MSA alignments and conservation scores."""
from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Optional


logger = logging.getLogger(__name__)


def save_alignment(alignment_fasta: Path, output_dir: Path, uniprot_id: str) -> Path:
    """Save alignment FASTA to standard location.

    Args:
        alignment_fasta: Source alignment FASTA file
        output_dir: Output directory (e.g., Data/processed/msa/)
        uniprot_id: UniProt accession (e.g., P04637)

    Returns:
        Path to saved alignment file
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"{uniprot_id}_alignment.fasta"

    # Copy alignment
    output_path.write_bytes(alignment_fasta.read_bytes())

    logger.info(f"Saved alignment to {output_path}")
    return output_path


def load_alignment(msa_dir: Path, uniprot_id: str) -> Optional[Path]:
    """Load cached alignment FASTA.

    Args:
        msa_dir: MSA directory (e.g., Data/processed/msa/)
        uniprot_id: UniProt accession

    Returns:
        Path to alignment file, or None if not found
    """
    alignment_path = msa_dir / f"{uniprot_id}_alignment.fasta"

    if not alignment_path.exists():
        logger.warning(f"Alignment not found: {alignment_path}")
        return None

    logger.info(f"Loaded alignment from {alignment_path}")
    return alignment_path


def save_conservation(
    conservation_scores: dict[int, float],
    output_dir: Path,
    uniprot_id: str,
) -> Path:
    """Save conservation scores to JSON.

    Args:
        conservation_scores: Position (1-indexed) → conservation score (0-1)
        output_dir: Output directory
        uniprot_id: UniProt accession

    Returns:
        Path to saved JSON file
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"{uniprot_id}_conservation.json"

    # Convert int keys to str for JSON
    json_data = {str(k): v for k, v in conservation_scores.items()}

    # Write JSON
    output_path.write_text(json.dumps(json_data, indent=2), encoding="utf-8")

    logger.info(f"Saved conservation scores to {output_path}")
    logger.info(f"  Positions: {len(conservation_scores)}")

    return output_path


def load_conservation(msa_dir: Path, uniprot_id: str) -> Optional[dict[int, float]]:
    """Load cached conservation scores.

    Args:
        msa_dir: MSA directory (e.g., Data/processed/msa/)
        uniprot_id: UniProt accession

    Returns:
        Position → conservation score, or None if not found
    """
    conservation_path = msa_dir / f"{uniprot_id}_conservation.json"

    if not conservation_path.exists():
        logger.warning(f"Conservation scores not found: {conservation_path}")
        return None

    # Load JSON
    json_data = json.loads(conservation_path.read_text(encoding="utf-8"))

    # Convert str keys back to int
    conservation_scores = {int(k): v for k, v in json_data.items()}

    logger.info(f"Loaded conservation scores from {conservation_path}")
    logger.info(f"  Positions: {len(conservation_scores)}")

    return conservation_scores


def get_or_create_conservation(
    msa_dir: Path,
    uniprot_id: str,
    query_seq: str,
    recompute: bool = False,
    **mafft_kwargs,
) -> dict[int, float]:
    """Get cached conservation scores or compute if missing.

    Args:
        msa_dir: MSA directory
        uniprot_id: UniProt accession
        query_seq: Query protein sequence
        recompute: Force recomputation even if cache exists
        **mafft_kwargs: Additional arguments for MAFFT (min_identity, max_seqs, etc.)

    Returns:
        Position → conservation score
    """
    # Check cache first
    if not recompute:
        conservation_scores = load_conservation(msa_dir, uniprot_id)
        if conservation_scores is not None:
            return conservation_scores

    logger.info("Computing conservation scores (cache miss or recompute=True)")

    # Import here to avoid circular dependency
    from .mafft_runner import fetch_homologs, run_mafft, validate_alignment
    from .conservation_scores import compute_conservation_scores

    import tempfile

    # Fetch homologs
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        homologs_fasta = tmpdir / "homologs.fasta"
        alignment_fasta = tmpdir / "alignment.fasta"

        # Fetch
        n_seqs = fetch_homologs(
            query_seq,
            uniprot_id,
            homologs_fasta,
            min_identity=mafft_kwargs.get("min_identity", 0.5),
            max_seqs=mafft_kwargs.get("max_seqs", 100),
            database=mafft_kwargs.get("database", "uniref90"),
        )

        if n_seqs < 2:
            logger.warning("Insufficient homologs, returning empty conservation scores")
            return {}

        # Align
        success = run_mafft(
            homologs_fasta,
            alignment_fasta,
            algorithm=mafft_kwargs.get("algorithm", "auto"),
            threads=mafft_kwargs.get("threads", 4),
        )

        if not success:
            logger.error("MAFFT alignment failed")
            return {}

        # Validate
        validation = validate_alignment(alignment_fasta)
        if not validation["valid"]:
            logger.error(f"Alignment validation failed: {validation['reason']}")
            return {}

        logger.info(f"Alignment quality:")
        logger.info(f"  Sequences: {validation['n_sequences']}")
        logger.info(f"  Length: {validation['alignment_length']}")
        logger.info(f"  Mean gap fraction: {validation['mean_gap_fraction']:.3f}")

        # Compute conservation
        conservation_scores = compute_conservation_scores(
            alignment_fasta,
            uniprot_id,
            max_gap_fraction=mafft_kwargs.get("max_gap_fraction", 0.3),
        )

        # Save cache
        save_alignment(alignment_fasta, msa_dir, uniprot_id)
        save_conservation(conservation_scores, msa_dir, uniprot_id)

    return conservation_scores
