"""MAFFT multiple sequence alignment runner with UniRef90 homolog fetching."""
from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path
from typing import Optional
import logging

from Bio import SeqIO, Entrez
from Bio.Blast import NCBIWWW, NCBIXML


logger = logging.getLogger(__name__)


def fetch_homologs(
    query_seq: str,
    uniprot_id: str,
    output_fasta: Path,
    min_identity: float = 0.5,
    max_seqs: int = 100,
    database: str = "uniref90",
) -> int:
    """Fetch homologous sequences from UniRef using BLAST.

    Args:
        query_seq: Query protein sequence (single-letter amino acids)
        uniprot_id: UniProt accession (e.g., P04637 for TP53)
        output_fasta: Path to save FASTA file with query + homologs
        min_identity: Minimum sequence identity threshold (0-1)
        min_identity: Minimum sequence identity threshold (0-1)
        max_seqs: Maximum number of homologs to fetch
        database: BLAST database (uniref90, uniref50, nr)

    Returns:
        Number of sequences written (including query)
    """
    logger.info(f"Fetching homologs for {uniprot_id} from {database}")
    logger.info(f"  Min identity: {min_identity:.1%}, Max seqs: {max_seqs}")

    # Write query sequence
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp:
        tmp_query = Path(tmp.name)
        tmp.write(f">{uniprot_id}\n{query_seq}\n")

    try:
        # Run BLAST search via web API
        logger.info("Running BLAST search (this may take 5-10 minutes)...")
        result_handle = NCBIWWW.qblast(
            program="blastp",
            database=database,
            sequence=query_seq,
            hitlist_size=max_seqs * 2,  # Fetch more to filter by identity
            expect=1e-10,
        )

        # Parse BLAST results
        blast_records = NCBIXML.parse(result_handle)
        record = next(blast_records)
        result_handle.close()

        # Filter alignments by identity
        homologs = []
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                identity = hsp.identities / hsp.align_length
                if identity >= min_identity:
                    # Extract accession from hit title
                    hit_id = alignment.title.split('|')[1] if '|' in alignment.title else alignment.hit_id
                    homologs.append({
                        'id': hit_id,
                        'identity': identity,
                        'evalue': hsp.expect,
                        'seq': hsp.sbjct.replace('-', ''),  # Remove gaps
                    })
                    break  # Take first HSP per alignment

        # Sort by identity and take top N
        homologs.sort(key=lambda x: x['identity'], reverse=True)
        homologs = homologs[:max_seqs]

        logger.info(f"Found {len(homologs)} homologs (identity ≥ {min_identity:.1%})")

        # Write FASTA file (query + homologs)
        with output_fasta.open('w') as f:
            # Write query first
            f.write(f">{uniprot_id}\n{query_seq}\n")

            # Write homologs
            for h in homologs:
                f.write(f">{h['id']} identity={h['identity']:.3f}\n")
                f.write(f"{h['seq']}\n")

        return len(homologs) + 1  # +1 for query

    except Exception as e:
        logger.error(f"BLAST search failed: {e}")
        # Fallback: write query-only FASTA
        with output_fasta.open('w') as f:
            f.write(f">{uniprot_id}\n{query_seq}\n")
        logger.warning("Created query-only FASTA (no homologs)")
        return 1

    finally:
        tmp_query.unlink(missing_ok=True)


def run_mafft(
    input_fasta: Path,
    output_fasta: Path,
    algorithm: str = "auto",
    threads: int = 4,
    timeout: int = 600,
) -> bool:
    """Run MAFFT multiple sequence alignment.

    Args:
        input_fasta: Path to input FASTA file (unaligned sequences)
        output_fasta: Path to output FASTA file (aligned sequences)
        algorithm: MAFFT algorithm (auto, linsi, ginsi, einsi, fftnsi, fftns)
        threads: Number of CPU threads to use
        timeout: Maximum runtime in seconds

    Returns:
        True if successful, False otherwise
    """
    logger.info(f"Running MAFFT alignment (algorithm={algorithm}, threads={threads})")

    # Build MAFFT command
    cmd = ["mafft"]

    # Algorithm selection
    if algorithm == "auto":
        cmd.append("--auto")
    elif algorithm == "linsi":
        cmd.append("--localpair")
        cmd.append("--maxiterate")
        cmd.append("1000")
    elif algorithm in ["ginsi", "einsi", "fftnsi", "fftns"]:
        cmd.append(f"--{algorithm}")
    else:
        logger.warning(f"Unknown algorithm '{algorithm}', using --auto")
        cmd.append("--auto")

    # Threading
    cmd.extend(["--thread", str(threads)])

    # Quiet mode (suppress progress messages)
    cmd.append("--quiet")

    # Input file
    cmd.append(str(input_fasta))

    try:
        # Run MAFFT
        logger.debug(f"Command: {' '.join(cmd)}")
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=timeout,
            text=True,
            check=True,
        )

        # Write output
        output_fasta.parent.mkdir(parents=True, exist_ok=True)
        output_fasta.write_text(result.stdout)

        # Count aligned sequences
        n_seqs = len(list(SeqIO.parse(output_fasta, "fasta")))
        logger.info(f"Alignment complete: {n_seqs} sequences")

        return True

    except subprocess.TimeoutExpired:
        logger.error(f"MAFFT timed out after {timeout} seconds")
        return False

    except subprocess.CalledProcessError as e:
        logger.error(f"MAFFT failed with exit code {e.returncode}")
        logger.error(f"stderr: {e.stderr}")
        return False

    except Exception as e:
        logger.error(f"Unexpected error running MAFFT: {e}")
        return False


def validate_alignment(alignment_fasta: Path, max_gap_fraction: float = 0.3) -> dict[str, any]:
    """Validate alignment quality.

    Args:
        alignment_fasta: Path to aligned FASTA file
        max_gap_fraction: Maximum allowed gap fraction per position

    Returns:
        Dictionary with validation metrics
    """
    sequences = list(SeqIO.parse(alignment_fasta, "fasta"))

    if not sequences:
        return {"valid": False, "reason": "No sequences found"}

    n_seqs = len(sequences)
    aln_length = len(sequences[0].seq)

    # Check all sequences have same length
    if not all(len(seq.seq) == aln_length for seq in sequences):
        return {"valid": False, "reason": "Sequences have different lengths"}

    # Calculate gap fraction per position
    gap_counts = [0] * aln_length
    for seq in sequences:
        for i, aa in enumerate(seq.seq):
            if aa == '-':
                gap_counts[i] += 1

    gap_fractions = [count / n_seqs for count in gap_counts]
    high_gap_positions = sum(1 for frac in gap_fractions if frac > max_gap_fraction)

    return {
        "valid": True,
        "n_sequences": n_seqs,
        "alignment_length": aln_length,
        "mean_gap_fraction": sum(gap_fractions) / len(gap_fractions),
        "max_gap_fraction": max(gap_fractions),
        "high_gap_positions": high_gap_positions,
        "high_gap_ratio": high_gap_positions / aln_length,
    }
