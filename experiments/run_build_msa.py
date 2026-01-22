#!/usr/bin/env python3
"""Build MSA and compute conservation scores for TP53.

Usage:
    python -m experiments.run_build_msa [--recompute]
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

# Add repo root to path
repo_root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(repo_root))

from src.utils.config_loader import load_config
from src.utils.logging_config import setup_logging
from src.msa.alignment_cache import get_or_create_conservation, load_conservation
from src.msa.conservation_scores import identify_conserved_regions


def _paths_from_config(paths_cfg: dict) -> dict[str, Path]:
    """Parse paths from config."""
    root = Path(paths_cfg.get("project_root", "."))
    data_root = root / paths_cfg.get("data_root", "Data")

    return {
        "root": root,
        "data": data_root,
        "processed": data_root / paths_cfg.get("data", {}).get("processed", "processed"),
        "msa": data_root / paths_cfg.get("data", {}).get("processed", "processed") / "msa",
    }


def main(recompute: bool = False):
    """Build MSA and conservation scores."""
    # Setup logging
    setup_logging()
    logger = logging.getLogger(__name__)

    logger.info("=" * 70)
    logger.info("MSA Builder: TP53 Conservation Analysis")
    logger.info("=" * 70)

    # Load config
    configs = load_config()
    paths = _paths_from_config(configs.get("paths", {}))
    p53_cfg = configs.get("p53", {})
    msa_cfg = p53_cfg.get("msa", {})

    # Get TP53 sequence and UniProt ID
    uniprot_id = p53_cfg.get("uniprot_id", "P04637")
    sequence = p53_cfg.get("sequence", "")

    # Fetch from UniProt if not in config
    if not sequence:
        logger.info(f"Sequence not in config, fetching from UniProt ({uniprot_id})...")
        try:
            from Bio import Entrez, SeqIO
            import urllib.request

            url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
            with urllib.request.urlopen(url) as response:
                fasta_text = response.read().decode('utf-8')

            from io import StringIO
            fasta_io = StringIO(fasta_text)
            record = next(SeqIO.parse(fasta_io, "fasta"))
            sequence = str(record.seq)

            logger.info(f"  Fetched sequence: {len(sequence)} residues")

        except Exception as e:
            logger.error(f"Failed to fetch sequence from UniProt: {e}")
            logger.error("Please add 'sequence' field to configs/p53.yaml")
            return 1

    logger.info(f"\nProtein: TP53 ({uniprot_id})")
    logger.info(f"Sequence length: {len(sequence)} residues")
    logger.info(f"Output directory: {paths['msa']}")

    # MSA parameters
    min_identity = msa_cfg.get("min_sequence_identity", 0.5)
    max_gap_fraction = msa_cfg.get("max_gap_fraction", 0.3)
    database = msa_cfg.get("database", "uniref90")

    logger.info(f"\nMSA Parameters:")
    logger.info(f"  Database: {database}")
    logger.info(f"  Min sequence identity: {min_identity:.1%}")
    logger.info(f"  Max gap fraction: {max_gap_fraction:.1%}")

    if recompute:
        logger.info(f"  Recompute: True (ignoring cache)")
    else:
        logger.info(f"  Recompute: False (using cache if available)")

    # Check if already computed
    if not recompute:
        existing = load_conservation(paths['msa'], uniprot_id)
        if existing is not None:
            logger.info(f"\n✓ Conservation scores already exist ({len(existing)} positions)")
            logger.info(f"  Use --recompute to regenerate")
            return 0

    # Compute conservation scores
    logger.info("\nStep 1: Fetching homologous sequences...")
    logger.info("  (This may take 5-10 minutes for BLAST search)")

    try:
        conservation_scores = get_or_create_conservation(
            msa_dir=paths['msa'],
            uniprot_id=uniprot_id,
            query_seq=sequence,
            recompute=recompute,
            min_identity=min_identity,
            max_gap_fraction=max_gap_fraction,
            database=database,
            max_seqs=100,
            algorithm="auto",
            threads=4,
        )

        if not conservation_scores:
            logger.error("Failed to compute conservation scores")
            return 1

        logger.info(f"\n{'=' * 70}")
        logger.info("✓ MSA and Conservation Scores Computed")
        logger.info(f"{'=' * 70}")

        # Statistics
        scores = list(conservation_scores.values())
        logger.info(f"\nConservation Statistics:")
        logger.info(f"  Positions analyzed: {len(scores)}")
        logger.info(f"  Mean conservation: {sum(scores) / len(scores):.3f}")
        logger.info(f"  Min conservation: {min(scores):.3f}")
        logger.info(f"  Max conservation: {max(scores):.3f}")

        # Conservation bins
        high = sum(1 for s in scores if s > 0.8)
        med = sum(1 for s in scores if 0.5 <= s <= 0.8)
        low = sum(1 for s in scores if s < 0.5)

        logger.info(f"\nConservation Distribution:")
        logger.info(f"  High (>0.8): {high} positions ({100*high/len(scores):.1f}%)")
        logger.info(f"  Medium (0.5-0.8): {med} positions ({100*med/len(scores):.1f}%)")
        logger.info(f"  Low (<0.5): {low} positions ({100*low/len(scores):.1f}%)")

        # Identify conserved regions
        regions = identify_conserved_regions(conservation_scores, threshold=0.8, min_length=3)

        logger.info(f"\nConserved Regions (>0.8, ≥3 residues):")
        if regions:
            for i, (start, end) in enumerate(regions, 1):
                length = end - start + 1
                mean_cons = sum(conservation_scores[pos] for pos in range(start, end+1)) / length
                logger.info(f"  {i}. Positions {start}-{end} (length={length}, mean={mean_cons:.3f})")
        else:
            logger.info("  No conserved regions found")

        # Output files
        logger.info(f"\nOutput Files:")
        alignment_path = paths['msa'] / f"{uniprot_id}_alignment.fasta"
        conservation_path = paths['msa'] / f"{uniprot_id}_conservation.json"

        if alignment_path.exists():
            logger.info(f"  Alignment: {alignment_path}")
        if conservation_path.exists():
            logger.info(f"  Conservation scores: {conservation_path}")

        logger.info(f"\nNext Steps:")
        logger.info(f"  1. Update configs/p53.yaml with:")
        logger.info(f"     msa:")
        logger.info(f"       enabled: true")
        logger.info(f"       precomputed_path: {conservation_path}")
        logger.info(f"  2. Run rescue design with MSA-based filtering")

        return 0

    except Exception as e:
        logger.error(f"Error computing conservation scores: {e}", exc_info=True)
        return 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build MSA and compute conservation scores for TP53")
    parser.add_argument(
        "--recompute",
        action="store_true",
        help="Force recomputation even if cache exists",
    )

    args = parser.parse_args()

    sys.exit(main(recompute=args.recompute))
