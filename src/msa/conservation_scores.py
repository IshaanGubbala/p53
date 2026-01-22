"""Shannon entropy-based conservation score calculation from MSA."""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional
import numpy as np

from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment


logger = logging.getLogger(__name__)


# Standard 20 amino acids
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


def shannon_entropy(column: list[str], pseudocount: float = 0.5) -> float:
    """Calculate Shannon entropy for an alignment column.

    Args:
        column: List of amino acids at a position (may include gaps '-')
        pseudocount: Pseudocount to add to frequencies (avoids log(0))

    Returns:
        Shannon entropy (0 = fully conserved, higher = more variable)
    """
    # Count amino acid frequencies (excluding gaps)
    aa_counts = {}
    n_non_gap = 0

    for aa in column:
        if aa != '-':
            aa_counts[aa] = aa_counts.get(aa, 0) + 1
            n_non_gap += 1

    if n_non_gap == 0:
        return 0.0  # All gaps

    # Add pseudocount and normalize
    freqs = []
    for aa in AMINO_ACIDS:
        count = aa_counts.get(aa, 0) + pseudocount
        freqs.append(count)

    total = sum(freqs)
    freqs = [f / total for f in freqs]

    # Calculate Shannon entropy: H = -sum(p * log2(p))
    entropy = -sum(p * np.log2(p) for p in freqs if p > 0)

    return entropy


def normalize_entropy(entropy: float, max_entropy: float = 4.32) -> float:
    """Normalize entropy to 0-1 scale (0 = conserved, 1 = variable).

    Args:
        entropy: Raw Shannon entropy
        max_entropy: Maximum possible entropy (log2(20) = 4.32 for 20 amino acids)

    Returns:
        Normalized entropy (0-1 scale)
    """
    return entropy / max_entropy


def compute_conservation_scores(
    alignment_fasta: Path,
    query_id: str,
    max_gap_fraction: float = 0.3,
) -> dict[int, float]:
    """Compute per-position conservation scores from MSA.

    Args:
        alignment_fasta: Path to aligned FASTA file
        query_id: ID of query sequence (used for residue numbering)
        max_gap_fraction: Positions with >this gap fraction get score=0

    Returns:
        Dictionary mapping position (1-indexed) to conservation score (0-1)
        Score: 0 = variable (low conservation), 1 = conserved (high conservation)
    """
    logger.info(f"Computing conservation scores from {alignment_fasta}")

    # Load alignment
    sequences = list(SeqIO.parse(alignment_fasta, "fasta"))

    if not sequences:
        logger.error("No sequences found in alignment")
        return {}

    # Find query sequence
    query_seq = None
    for seq in sequences:
        if query_id in seq.id:
            query_seq = seq
            break

    if query_seq is None:
        logger.error(f"Query sequence '{query_id}' not found in alignment")
        return {}

    n_seqs = len(sequences)
    aln_length = len(query_seq.seq)

    logger.info(f"Alignment: {n_seqs} sequences, {aln_length} positions")

    # Extract columns
    columns = []
    for i in range(aln_length):
        column = [str(seq.seq[i]).upper() for seq in sequences]
        columns.append(column)

    # Compute conservation scores
    conservation_scores = {}
    query_pos = 0  # Current position in query (ungapped)

    for i, column in enumerate(columns):
        # Skip if query has gap at this position
        if str(query_seq.seq[i]) == '-':
            continue

        query_pos += 1  # Increment ungapped position

        # Calculate gap fraction
        gap_count = sum(1 for aa in column if aa == '-')
        gap_frac = gap_count / n_seqs

        # High-gap positions get score = 0 (uncertain)
        if gap_frac > max_gap_fraction:
            conservation_scores[query_pos] = 0.0
            continue

        # Calculate Shannon entropy
        entropy = shannon_entropy(column)

        # Normalize to 0-1 scale
        norm_entropy = normalize_entropy(entropy)

        # Invert so 1 = conserved, 0 = variable
        conservation = 1.0 - norm_entropy

        conservation_scores[query_pos] = conservation

    logger.info(f"Computed conservation for {len(conservation_scores)} positions")
    logger.info(f"  Mean conservation: {np.mean(list(conservation_scores.values())):.3f}")
    logger.info(f"  Highly conserved (>0.8): {sum(1 for v in conservation_scores.values() if v > 0.8)}")

    return conservation_scores


def identify_conserved_regions(
    conservation_scores: dict[int, float],
    threshold: float = 0.8,
    min_length: int = 3,
) -> list[tuple[int, int]]:
    """Identify contiguous conserved regions.

    Args:
        conservation_scores: Position → conservation score (0-1)
        threshold: Minimum conservation score to consider
        min_length: Minimum region length

    Returns:
        List of (start, end) tuples for conserved regions (inclusive, 1-indexed)
    """
    positions = sorted(conservation_scores.keys())

    regions = []
    region_start = None

    for pos in positions:
        score = conservation_scores[pos]

        if score >= threshold:
            if region_start is None:
                region_start = pos
        else:
            if region_start is not None:
                region_end = positions[positions.index(pos) - 1]
                if (region_end - region_start + 1) >= min_length:
                    regions.append((region_start, region_end))
                region_start = None

    # Handle last region
    if region_start is not None:
        region_end = positions[-1]
        if (region_end - region_start + 1) >= min_length:
            regions.append((region_start, region_end))

    logger.info(f"Identified {len(regions)} conserved regions (threshold={threshold})")

    return regions
