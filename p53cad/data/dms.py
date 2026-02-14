from __future__ import annotations

import pandas as pd
import numpy as np
import re
from pathlib import Path
from typing import Optional, Tuple, List, Dict
from p53cad.core.logging import get_logger

logger = get_logger(__name__)

# =============================================================================
# DATA SOURCES
# =============================================================================
# Full Giacomelli 2018 saturation mutagenesis data (~8260 variants)
# Paper: https://www.nature.com/articles/s41588-018-0204-y
#
# Score columns available:
# - A549_p53WT_Nutlin-3_Z-score: Main functional score (p53 activation assay)
# - A549_p53NULL_Nutlin-3_Z-score: Control (no p53)
# - A549_p53NULL_Etoposide_Z-score: DNA damage response
# =============================================================================

# Wild-type p53 sequence (393 amino acids)
P53_WT = (
    "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP"
    "DEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAK"
    "SVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHE"
    "RCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNS"
    "SCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELP"
    "PGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPG"
    "GSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD"
)

# Default path to full saturation DMS data
DEFAULT_DMS_PATH = Path(__file__).parent.parent.parent / "data" / "raw" / "p53_DMS_Giacomelli_2018_FULL.csv"
# Fallback to old cell-line data if full data not available
LEGACY_DMS_PATH = Path(__file__).parent.parent.parent / "data" / "raw" / "p53_DMS_Giacomelli_2018.csv"


def parse_single_mutation(mutation_str: str) -> Optional[Tuple[str, int, str]]:
    """
    Parse a single mutation string like 'R175H' into (wt_aa, position, mut_aa).

    Returns None for complex mutations (frameshift, nonsense, splice, deletions).
    """
    mutation_str = mutation_str.strip()

    # Skip complex mutations
    if any(x in mutation_str for x in ['fs', '*', 'del', 'ins', 'splice', '>']):
        return None

    # Match standard format: A123B (single letter, number, single letter)
    match = re.match(r'^([A-Z])(\d+)([A-Z])$', mutation_str)
    if match:
        wt_aa, pos_str, mut_aa = match.groups()
        pos = int(pos_str)
        # Validate position is within p53 sequence
        if 1 <= pos <= len(P53_WT):
            return (wt_aa, pos, mut_aa)

    return None


def parse_mutation(mutation: str) -> List[Tuple[str, int, str]]:
    """
    Parse mutation string which may contain multiple mutations (comma-separated).

    Returns list of (wt_aa, position, mut_aa) tuples for valid point mutations.
    Returns empty list if any mutation is invalid/complex.
    """
    # Handle quoted strings from CSV
    mutation = mutation.strip('"').strip()

    # Split by comma for double/multiple mutants
    parts = [p.strip() for p in mutation.split(',')]

    parsed = []
    for part in parts:
        result = parse_single_mutation(part)
        if result is None:
            return []  # If any part is invalid, skip entire entry
        parsed.append(result)

    return parsed


def apply_mutation(wt_seq: str, mutation: str) -> Optional[str]:
    """
    Applies a mutation string (e.g. R175H or 'R175H, T211I') to WT sequence.

    Returns None if mutation is complex (fs, *, del) or invalid.
    """
    parsed = parse_mutation(mutation)
    if not parsed:
        return None

    seq = list(wt_seq)

    for wt_aa, pos, mut_aa in parsed:
        idx = pos - 1  # Convert 1-indexed to 0-indexed

        # Validate WT amino acid matches (with some tolerance for data errors)
        if idx < len(seq):
            # Note: Some mutations in the dataset may have WT mismatches due to
            # different reference sequences or data entry errors. We still apply them.
            seq[idx] = mut_aa
        else:
            return None

    return "".join(seq)


def load_dms_data(file_path: Optional[Path | str] = None, score_column: str = "A549_p53WT_Nutlin-3_Z-score") -> pd.DataFrame:
    """
    Load Giacomelli 2018 p53 DMS data.

    Supports two formats:
    1. Full saturation format (8260 entries): Allele, AA_wt, AA_variant, Position, Z-scores
    2. Legacy cell-line format (529 entries): mutation, score

    The main functional score is A549_p53WT_Nutlin-3_Z-score which measures
    p53 transcriptional activity. Negative = loss of function, positive = functional.

    Args:
        file_path: Path to DMS CSV. Uses default path if not provided.
        score_column: Which Z-score column to use (default: A549_p53WT_Nutlin-3_Z-score)

    Returns:
        DataFrame with columns: mutation, score, pos, wt, alt, sequence
    """
    path = Path(file_path) if file_path else DEFAULT_DMS_PATH

    # Try paths in order: specified > default full > legacy > parent dir
    if not path.exists():
        if LEGACY_DMS_PATH.exists():
            path = LEGACY_DMS_PATH
            logger.info("Using legacy DMS data (cell-line format)")
        else:
            alt_path = Path(__file__).parent.parent.parent.parent / "p53_DMS_Giacomelli_2018.csv"
            if alt_path.exists():
                path = alt_path
            else:
                raise FileNotFoundError(
                    f"DMS data not found. Searched: {DEFAULT_DMS_PATH}, {LEGACY_DMS_PATH}"
                )

    logger.info(f"Loading DMS data from {path}")
    df = pd.read_csv(path)

    # Detect format based on columns
    is_saturation_format = 'Allele' in df.columns and 'Position' in df.columns
    is_legacy_format = 'mutation' in df.columns and 'score' in df.columns

    if is_saturation_format:
        return _load_saturation_format(df, score_column)
    elif is_legacy_format:
        return _load_legacy_format(df)
    else:
        raise ValueError(f"Unknown DMS file format. Columns: {df.columns.tolist()}")


def _load_saturation_format(df: pd.DataFrame, score_column: str) -> pd.DataFrame:
    """Load full saturation mutagenesis format (8260 entries).

    NOTE: Nutlin-3 Z-scores are INVERTED - positive = loss of function.
    We negate scores so negative = pathogenic (matches clinical convention).
    """
    logger.info(f"Detected saturation format with {len(df)} entries")

    if score_column not in df.columns:
        # Try to find a suitable score column
        score_cols = [c for c in df.columns if 'Z-score' in c or 'score' in c.lower()]
        if score_cols:
            score_column = score_cols[0]
            logger.warning(f"Score column not found, using: {score_column}")
        else:
            raise ValueError(f"No score column found. Available: {df.columns.tolist()}")

    processed_data = []
    skipped = 0

    for _, row in df.iterrows():
        allele = str(row.get('Allele', ''))
        wt_aa = str(row.get('AA_wt', ''))
        mut_aa = str(row.get('AA_variant', ''))
        pos = int(row.get('Position', 0))
        raw_score = row.get(score_column)

        # Skip rows with missing data or stop codons (Z)
        if pd.isna(raw_score) or not allele or mut_aa == 'Z' or pos == 0:
            skipped += 1
            continue

        # NEGATE score: Nutlin-3 positive = loss of function, we want negative = pathogenic
        score = -float(raw_score)

        # Build mutation string (e.g., "R175H")
        mutation = f"{wt_aa}{pos}{mut_aa}"

        # Generate mutant sequence
        sequence = apply_mutation(P53_WT, mutation)
        if not sequence:
            skipped += 1
            continue

        processed_data.append({
            'mutation': mutation,
            'score': score,
            'pos': pos,
            'wt': wt_aa,
            'alt': mut_aa,
            'sequence': sequence,
            'n_mutations': 1,
            'is_hotspot': pos in [175, 248, 273, 245, 249, 282, 220],
            'allele': allele  # Keep original allele notation
        })

    result_df = pd.DataFrame(processed_data)
    logger.info(f"Loaded {len(result_df)} valid saturation DMS entries (skipped {skipped})")
    logger.info(f"  - Positions covered: {result_df['pos'].nunique()}")
    logger.info(f"  - Score range: [{result_df['score'].min():.2f}, {result_df['score'].max():.2f}]")

    return result_df


def _load_legacy_format(df: pd.DataFrame) -> pd.DataFrame:
    """Load legacy cell-line format (529 entries)."""
    logger.info(f"Detected legacy format with {len(df)} entries")

    processed_data = []
    skipped_complex = 0
    skipped_invalid = 0

    for _, row in df.iterrows():
        mutation = str(row['mutation'])
        score = float(row['score'])

        parsed = parse_mutation(mutation)
        if not parsed:
            skipped_complex += 1
            continue

        if len(parsed) == 1:
            wt_aa, pos, mut_aa = parsed[0]
            sequence = apply_mutation(P53_WT, mutation)

            if sequence:
                processed_data.append({
                    'mutation': mutation,
                    'score': score,
                    'pos': pos,
                    'wt': wt_aa,
                    'alt': mut_aa,
                    'sequence': sequence,
                    'n_mutations': 1,
                    'is_hotspot': pos in [175, 248, 273, 245, 249, 282, 220]
                })
            else:
                skipped_invalid += 1
        else:
            sequence = apply_mutation(P53_WT, mutation)
            if sequence:
                wt_aa, pos, mut_aa = parsed[0]
                processed_data.append({
                    'mutation': mutation,
                    'score': score,
                    'pos': pos,
                    'wt': wt_aa,
                    'alt': mut_aa,
                    'sequence': sequence,
                    'n_mutations': len(parsed),
                    'is_hotspot': any(p[1] in [175, 248, 273, 245, 249, 282, 220] for p in parsed)
                })

    result_df = pd.DataFrame(processed_data)
    logger.info(f"Loaded {len(result_df)} valid legacy DMS entries")
    logger.info(f"  - Single mutations: {len(result_df[result_df['n_mutations'] == 1])}")
    logger.info(f"  - Multiple mutations: {len(result_df[result_df['n_mutations'] > 1])}")
    logger.info(f"  - Skipped complex (fs/*/del): {skipped_complex}")
    logger.info(f"  - Skipped invalid: {skipped_invalid}")
    logger.info(f"  - Hotspot mutations: {len(result_df[result_df['is_hotspot']])}")

    return result_df


def hydrate_sequences(df: pd.DataFrame) -> pd.DataFrame:
    """
    Adds a 'sequence' column to a dataframe with a 'mutation' column.
    For pre-processed data that already has sequences, returns as-is.
    """
    if 'sequence' in df.columns:
        return df

    hydrated = []
    for _, row in df.iterrows():
        seq = apply_mutation(P53_WT, str(row["mutation"]))
        if seq:
            new_row = row.to_dict()
            new_row["sequence"] = seq
            hydrated.append(new_row)

    return pd.DataFrame(hydrated)


def get_mutation_statistics(df: pd.DataFrame) -> dict:
    """
    Compute statistics about the DMS dataset.
    """
    stats = {
        'total_entries': len(df),
        'unique_positions': df['pos'].nunique() if 'pos' in df.columns else 0,
        'score_mean': df['score'].mean(),
        'score_std': df['score'].std(),
        'score_min': df['score'].min(),
        'score_max': df['score'].max(),
        'pathogenic_count': len(df[df['score'] < -0.5]),
        'neutral_count': len(df[(df['score'] >= -0.5) & (df['score'] <= 0.5)]),
        'functional_count': len(df[df['score'] > 0.5]),
    }

    if 'is_hotspot' in df.columns:
        stats['hotspot_count'] = len(df[df['is_hotspot']])
        stats['hotspot_mean_score'] = df[df['is_hotspot']]['score'].mean()

    return stats


def load_double_mutant_data(path: Optional[Path] = None) -> pd.DataFrame:
    """Load double-mutant DMS data from external source.

    Expected CSV columns: pos1, aa1, pos2, aa2, score (or mutation column
    with comma-separated pairs like 'R175H,N268D').

    Falls back to ``data/raw/p53_DMS_double_mutants.csv`` when *path* is None.
    Returns an empty DataFrame if the file is not found (graceful degradation).
    """
    default_path = Path(__file__).parent.parent.parent / "data" / "raw" / "p53_DMS_double_mutants.csv"
    path = Path(path) if path else default_path
    if not path.exists():
        logger.info("Double-mutant DMS file not found at %s; using additive fallback", path)
        return pd.DataFrame()

    logger.info("Loading double-mutant DMS data from %s", path)
    df = pd.read_csv(path)

    # Normalise columns
    if "mutation" in df.columns and "pos1" not in df.columns:
        # Parse 'R175H,N268D' format
        rows = []
        score_col = next((c for c in df.columns if "score" in c.lower() or "z" in c.lower()), None)
        if score_col is None:
            logger.warning("No score column found in double-mutant CSV")
            return pd.DataFrame()
        for _, row in df.iterrows():
            parts = parse_mutation(str(row["mutation"]))
            if len(parts) == 2:
                (_, p1, a1), (_, p2, a2) = parts
                rows.append({"pos1": p1, "aa1": a1, "pos2": p2, "aa2": a2, "score": float(row[score_col])})
        df = pd.DataFrame(rows)

    required = {"pos1", "aa1", "pos2", "aa2", "score"}
    if not required.issubset(set(df.columns)):
        logger.warning("Double-mutant CSV missing columns: %s", required - set(df.columns))
        return pd.DataFrame()

    logger.info("Loaded %d double-mutant DMS entries", len(df))
    return df


def get_pairwise_epistasis_lookup(
    dms_df: Optional[pd.DataFrame] = None,
    double_df: Optional[pd.DataFrame] = None,
) -> Dict[tuple, float]:
    """Build (pos1, aa1, pos2, aa2) → epistasis score lookup.

    Uses real double-mutant data when available.  Falls back to an additive
    model: epistasis = observed_double - (single_A + single_B).  When no
    double-mutant data exists, returns an empty dict (callers should fall
    back to structural heuristics).

    Returns
    -------
    dict
        Keys are ``(pos1, aa1, pos2, aa2)`` tuples (positions 1-indexed).
        Values are Z-score differences (positive = synergistic).
    """
    if double_df is None:
        double_df = load_double_mutant_data()
    if double_df.empty:
        return {}

    # Load single-mutant data for additive baseline
    if dms_df is None:
        try:
            dms_df = get_dms_data()
        except FileNotFoundError:
            dms_df = pd.DataFrame()

    single_lookup: Dict[tuple, float] = {}
    if not dms_df.empty and "pos" in dms_df.columns and "alt" in dms_df.columns and "score" in dms_df.columns:
        for _, row in dms_df.iterrows():
            single_lookup[(int(row["pos"]), str(row["alt"]))] = float(row["score"])

    result: Dict[tuple, float] = {}
    for _, row in double_df.iterrows():
        p1, a1 = int(row["pos1"]), str(row["aa1"]).upper()
        p2, a2 = int(row["pos2"]), str(row["aa2"]).upper()
        observed = float(row["score"])

        # Additive expectation (sum of singles)
        s1 = single_lookup.get((p1, a1), 0.0)
        s2 = single_lookup.get((p2, a2), 0.0)
        epistasis = observed - (s1 + s2)

        # Store both orderings for easy lookup
        result[(p1, a1, p2, a2)] = epistasis
        result[(p2, a2, p1, a1)] = epistasis

    logger.info("Built pairwise epistasis lookup: %d pairs", len(result) // 2)
    return result


# Convenience function for quick loading
def get_dms_data() -> pd.DataFrame:
    """Load DMS data with default settings."""
    return load_dms_data()


# =============================================================================
# PHYSICS-BASED FALLBACK SCORING
# =============================================================================
# When DMS data is unavailable, use amino acid properties to estimate impact

# Kyte-Doolittle hydropathy scale
HYDROPATHY = {
    'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5, 'M': 1.9, 'A': 1.8,
    'G': -0.4, 'T': -0.7, 'S': -0.8, 'W': -0.9, 'Y': -1.3, 'P': -1.6,
    'H': -3.2, 'E': -3.5, 'Q': -3.5, 'D': -3.5, 'N': -3.5, 'K': -3.9, 'R': -4.5
}

# Amino acid molecular volumes (Å³)
AA_VOLUME = {
    'G': 60, 'A': 88, 'S': 89, 'C': 108, 'D': 111, 'P': 112, 'N': 114,
    'T': 116, 'E': 138, 'V': 140, 'Q': 143, 'H': 153, 'M': 162, 'I': 166,
    'L': 166, 'K': 168, 'R': 173, 'F': 189, 'Y': 193, 'W': 227
}

# Amino acid charges at physiological pH
AA_CHARGE = {
    'K': 1, 'R': 1, 'H': 0.1,  # Basic (positive)
    'D': -1, 'E': -1,          # Acidic (negative)
    'A': 0, 'C': 0, 'F': 0, 'G': 0, 'I': 0, 'L': 0, 'M': 0,
    'N': 0, 'P': 0, 'Q': 0, 'S': 0, 'T': 0, 'V': 0, 'W': 0, 'Y': 0
}

# BLOSUM62 substitution matrix (symmetric)
# Higher = more conservative substitution
BLOSUM62 = {
    ('A', 'A'): 4, ('A', 'R'): -1, ('A', 'N'): -2, ('A', 'D'): -2, ('A', 'C'): 0,
    ('A', 'Q'): -1, ('A', 'E'): -1, ('A', 'G'): 0, ('A', 'H'): -2, ('A', 'I'): -1,
    ('A', 'L'): -1, ('A', 'K'): -1, ('A', 'M'): -1, ('A', 'F'): -2, ('A', 'P'): -1,
    ('A', 'S'): 1, ('A', 'T'): 0, ('A', 'W'): -3, ('A', 'Y'): -2, ('A', 'V'): 0,
    ('R', 'R'): 5, ('R', 'N'): 0, ('R', 'D'): -2, ('R', 'C'): -3, ('R', 'Q'): 1,
    ('R', 'E'): 0, ('R', 'G'): -2, ('R', 'H'): 0, ('R', 'I'): -3, ('R', 'L'): -2,
    ('R', 'K'): 2, ('R', 'M'): -1, ('R', 'F'): -3, ('R', 'P'): -2, ('R', 'S'): -1,
    ('R', 'T'): -1, ('R', 'W'): -3, ('R', 'Y'): -2, ('R', 'V'): -3,
    ('N', 'N'): 6, ('N', 'D'): 1, ('N', 'C'): -3, ('N', 'Q'): 0, ('N', 'E'): 0,
    ('N', 'G'): 0, ('N', 'H'): 1, ('N', 'I'): -3, ('N', 'L'): -3, ('N', 'K'): 0,
    ('N', 'M'): -2, ('N', 'F'): -3, ('N', 'P'): -2, ('N', 'S'): 1, ('N', 'T'): 0,
    ('N', 'W'): -4, ('N', 'Y'): -2, ('N', 'V'): -3,
    ('D', 'D'): 6, ('D', 'C'): -3, ('D', 'Q'): 0, ('D', 'E'): 2, ('D', 'G'): -1,
    ('D', 'H'): -1, ('D', 'I'): -3, ('D', 'L'): -4, ('D', 'K'): -1, ('D', 'M'): -3,
    ('D', 'F'): -3, ('D', 'P'): -1, ('D', 'S'): 0, ('D', 'T'): -1, ('D', 'W'): -4,
    ('D', 'Y'): -3, ('D', 'V'): -3,
    ('C', 'C'): 9, ('C', 'Q'): -3, ('C', 'E'): -4, ('C', 'G'): -3, ('C', 'H'): -3,
    ('C', 'I'): -1, ('C', 'L'): -1, ('C', 'K'): -3, ('C', 'M'): -1, ('C', 'F'): -2,
    ('C', 'P'): -3, ('C', 'S'): -1, ('C', 'T'): -1, ('C', 'W'): -2, ('C', 'Y'): -2,
    ('C', 'V'): -1,
    ('Q', 'Q'): 5, ('Q', 'E'): 2, ('Q', 'G'): -2, ('Q', 'H'): 0, ('Q', 'I'): -3,
    ('Q', 'L'): -2, ('Q', 'K'): 1, ('Q', 'M'): 0, ('Q', 'F'): -3, ('Q', 'P'): -1,
    ('Q', 'S'): 0, ('Q', 'T'): -1, ('Q', 'W'): -2, ('Q', 'Y'): -1, ('Q', 'V'): -2,
    ('E', 'E'): 5, ('E', 'G'): -2, ('E', 'H'): 0, ('E', 'I'): -3, ('E', 'L'): -3,
    ('E', 'K'): 1, ('E', 'M'): -2, ('E', 'F'): -3, ('E', 'P'): -1, ('E', 'S'): 0,
    ('E', 'T'): -1, ('E', 'W'): -3, ('E', 'Y'): -2, ('E', 'V'): -2,
    ('G', 'G'): 6, ('G', 'H'): -2, ('G', 'I'): -4, ('G', 'L'): -4, ('G', 'K'): -2,
    ('G', 'M'): -3, ('G', 'F'): -3, ('G', 'P'): -2, ('G', 'S'): 0, ('G', 'T'): -2,
    ('G', 'W'): -2, ('G', 'Y'): -3, ('G', 'V'): -3,
    ('H', 'H'): 8, ('H', 'I'): -3, ('H', 'L'): -3, ('H', 'K'): -1, ('H', 'M'): -2,
    ('H', 'F'): -1, ('H', 'P'): -2, ('H', 'S'): -1, ('H', 'T'): -2, ('H', 'W'): -2,
    ('H', 'Y'): 2, ('H', 'V'): -3,
    ('I', 'I'): 4, ('I', 'L'): 2, ('I', 'K'): -3, ('I', 'M'): 1, ('I', 'F'): 0,
    ('I', 'P'): -3, ('I', 'S'): -2, ('I', 'T'): -1, ('I', 'W'): -3, ('I', 'Y'): -1,
    ('I', 'V'): 3,
    ('L', 'L'): 4, ('L', 'K'): -2, ('L', 'M'): 2, ('L', 'F'): 0, ('L', 'P'): -3,
    ('L', 'S'): -2, ('L', 'T'): -1, ('L', 'W'): -2, ('L', 'Y'): -1, ('L', 'V'): 1,
    ('K', 'K'): 5, ('K', 'M'): -1, ('K', 'F'): -3, ('K', 'P'): -1, ('K', 'S'): 0,
    ('K', 'T'): -1, ('K', 'W'): -3, ('K', 'Y'): -2, ('K', 'V'): -2,
    ('M', 'M'): 5, ('M', 'F'): 0, ('M', 'P'): -2, ('M', 'S'): -1, ('M', 'T'): -1,
    ('M', 'W'): -1, ('M', 'Y'): -1, ('M', 'V'): 1,
    ('F', 'F'): 6, ('F', 'P'): -4, ('F', 'S'): -2, ('F', 'T'): -2, ('F', 'W'): 1,
    ('F', 'Y'): 3, ('F', 'V'): -1,
    ('P', 'P'): 7, ('P', 'S'): -1, ('P', 'T'): -1, ('P', 'W'): -4, ('P', 'Y'): -3,
    ('P', 'V'): -2,
    ('S', 'S'): 4, ('S', 'T'): 1, ('S', 'W'): -3, ('S', 'Y'): -2, ('S', 'V'): -2,
    ('T', 'T'): 5, ('T', 'W'): -2, ('T', 'Y'): -2, ('T', 'V'): 0,
    ('W', 'W'): 11, ('W', 'Y'): 2, ('W', 'V'): -3,
    ('Y', 'Y'): 7, ('Y', 'V'): -1,
    ('V', 'V'): 4,
}

# p53 DNA-binding domain residues (critical for function)
P53_DBD_RANGE = (94, 292)  # Core DNA-binding domain

# Known critical functional positions in p53
P53_CRITICAL_POSITIONS = {
    # Zinc coordination sites (critical for structure)
    176, 179, 238, 242,
    # DNA-contacting residues
    248, 273, 280, 241, 282,
    # L2 loop (involved in DNA binding)
    163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175,
    # L3 loop (involved in DNA binding)
    237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250,
    # LSH motif
    119, 120, 121
}


def get_blosum62_score(aa1: str, aa2: str) -> int:
    """Get BLOSUM62 substitution score (symmetric lookup)."""
    if aa1 == aa2:
        return BLOSUM62.get((aa1, aa1), 0)
    # Try both orderings since matrix is symmetric
    return BLOSUM62.get((aa1, aa2), BLOSUM62.get((aa2, aa1), -4))


def physics_based_score(mutation: str) -> Dict[str, float]:
    """
    Calculate physics-based mutation impact score when DMS data unavailable.

    Uses amino acid properties to estimate functional impact:
    - Hydrophobicity change (burial disruption)
    - Size change (steric clashes)
    - Charge change (electrostatic effects)
    - BLOSUM62 evolutionary score (conservation)
    - Position in critical regions

    Returns dict with component scores and combined estimate.
    Score interpretation: negative = likely deleterious, positive = likely tolerated
    """
    parsed = parse_single_mutation(mutation)
    if not parsed:
        return {'score': 0.0, 'confidence': 0.0, 'method': 'invalid_mutation'}

    wt_aa, pos, mut_aa = parsed

    # Component calculations
    hydro_wt = HYDROPATHY.get(wt_aa, 0)
    hydro_mut = HYDROPATHY.get(mut_aa, 0)
    hydro_change = abs(hydro_mut - hydro_wt)  # Larger = more disruptive

    vol_wt = AA_VOLUME.get(wt_aa, 120)
    vol_mut = AA_VOLUME.get(mut_aa, 120)
    vol_change = abs(vol_mut - vol_wt)  # Larger = more steric clash risk

    charge_wt = AA_CHARGE.get(wt_aa, 0)
    charge_mut = AA_CHARGE.get(mut_aa, 0)
    charge_change = abs(charge_mut - charge_wt)  # Charge reversal is severe

    blosum = get_blosum62_score(wt_aa, mut_aa)  # Higher = more conservative

    # Position-based modifiers
    in_dbd = P53_DBD_RANGE[0] <= pos <= P53_DBD_RANGE[1]
    is_critical = pos in P53_CRITICAL_POSITIONS

    # Calculate combined score (calibrated to DMS scale: -2 to +2)
    # More negative = more likely deleterious

    # Start with BLOSUM as base (normalized from -4 to +11 range)
    base_score = (blosum + 4) / 15 * 2 - 1  # Maps to roughly -1 to +1

    # Penalize large property changes
    hydro_penalty = -0.15 * hydro_change  # Max ~1.35 penalty
    vol_penalty = -0.005 * vol_change      # Max ~0.8 penalty
    charge_penalty = -0.5 * charge_change  # Max 1.0 penalty for charge reversal

    # Position modifiers
    position_modifier = 0.0
    if is_critical:
        position_modifier = -0.5  # Critical positions more sensitive
    elif in_dbd:
        position_modifier = -0.2  # DBD generally more sensitive

    # Combined score
    combined = base_score + hydro_penalty + vol_penalty + charge_penalty + position_modifier
    combined = np.clip(combined, -2.0, 1.5)  # Clamp to reasonable DMS range

    # Confidence based on how many property changes occur
    n_changes = sum([
        hydro_change > 2.0,
        vol_change > 30,
        charge_change > 0,
        blosum < 0
    ])
    confidence = 0.3 + 0.15 * n_changes  # Higher confidence when clear signal

    return {
        'score': float(combined),
        'confidence': float(confidence),
        'method': 'physics_based',
        'components': {
            'blosum62': blosum,
            'hydropathy_change': float(hydro_change),
            'volume_change': float(vol_change),
            'charge_change': float(charge_change),
            'in_dbd': in_dbd,
            'is_critical': is_critical,
            'base_score': float(base_score),
            'hydro_penalty': float(hydro_penalty),
            'vol_penalty': float(vol_penalty),
            'charge_penalty': float(charge_penalty),
            'position_modifier': float(position_modifier)
        }
    }


def get_mutation_score(mutation: str, dms_df: Optional[pd.DataFrame] = None) -> Dict[str, float]:
    """
    Get mutation score, preferring DMS data but falling back to physics-based estimate.

    Args:
        mutation: Mutation string like 'R175H'
        dms_df: Optional pre-loaded DMS dataframe

    Returns:
        Dict with score, confidence, and method used
    """
    # Try DMS lookup first
    if dms_df is not None and not dms_df.empty:
        match = dms_df[dms_df['mutation'] == mutation]
        if not match.empty:
            return {
                'score': float(match['score'].iloc[0]),
                'confidence': 1.0,  # High confidence for experimental data
                'method': 'dms_experimental'
            }

    # Fall back to physics-based score
    return physics_based_score(mutation)


def get_coverage_report(mutations: List[str], dms_df: Optional[pd.DataFrame] = None) -> Dict:
    """
    Generate a detailed coverage report showing DMS vs physics-based scoring.

    Useful for understanding what data IS available vs estimated.
    """
    if dms_df is None:
        try:
            dms_df = get_dms_data()
        except FileNotFoundError:
            dms_df = pd.DataFrame()

    dms_mutations = set(dms_df['mutation'].tolist()) if not dms_df.empty else set()
    dms_positions = set(dms_df['pos'].tolist()) if 'pos' in dms_df.columns else set()

    results = []
    dms_count = 0
    physics_count = 0

    for mut in mutations:
        parsed = parse_single_mutation(mut)
        if not parsed:
            continue

        wt_aa, pos, mut_aa = parsed

        if mut in dms_mutations:
            score_info = {'score': float(dms_df[dms_df['mutation'] == mut]['score'].iloc[0]),
                         'confidence': 1.0, 'method': 'dms_experimental'}
            dms_count += 1
        else:
            score_info = physics_based_score(mut)
            physics_count += 1

        results.append({
            'mutation': mut,
            'position': pos,
            'wt_aa': wt_aa,
            'mut_aa': mut_aa,
            'score': score_info['score'],
            'confidence': score_info['confidence'],
            'method': score_info['method'],
            'position_has_dms': pos in dms_positions
        })

    return {
        'mutations': results,
        'summary': {
            'total': len(results),
            'dms_experimental': dms_count,
            'physics_estimated': physics_count,
            'dms_coverage': dms_count / max(len(results), 1),
            'positions_with_any_dms': sum(1 for r in results if r['position_has_dms']),
            'dms_dataset_size': len(dms_df),
            'dms_positions_covered': len(dms_positions)
        },
        'data_sources': {
            'current_file': 'Giacomelli 2018 cell line data (~367 entries)',
            'full_saturation': 'Download from MaveDB: https://www.mavedb.org/',
            'alternative': 'Nature Genetics 2024: https://www.nature.com/articles/s41588-024-02039-4'
        }
    }


def print_data_source_info() -> None:
    """Print helpful information about getting more DMS data."""
    info = """
╔══════════════════════════════════════════════════════════════════════╗
║                    p53 DMS DATA SOURCES                              ║
╠══════════════════════════════════════════════════════════════════════╣
║ CURRENT DATA: Giacomelli 2018 cell line (~367 entries, 88 positions) ║
║                                                                      ║
║ FOR FULL SATURATION MUTAGENESIS DATA (~8000+ variants):              ║
║                                                                      ║
║ 1. MaveDB (Recommended)                                              ║
║    → https://www.mavedb.org/                                         ║
║    → Search "TP53" or "p53 Giacomelli"                               ║
║    → Download CSV with all variant effects                           ║
║                                                                      ║
║ 2. Nature Genetics 2024 - Deep CRISPR Study (Most Comprehensive)     ║
║    → https://www.nature.com/articles/s41588-024-02039-4              ║
║    → Supplementary data contains all TP53 functional data            ║
║                                                                      ║
║ 3. Original Giacomelli 2018 Paper                                    ║
║    → https://www.nature.com/articles/s41588-018-0204-y               ║
║    → Supplementary Tables in paper                                   ║
║                                                                      ║
║ Place downloaded CSV in: data/raw/p53_DMS_Giacomelli_2018.csv        ║
╚══════════════════════════════════════════════════════════════════════╝
"""
    print(info)


def get_position_statistics(dms_df: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """
    Get per-position statistics from available DMS data.

    Shows what positions HAVE data and their score distributions.
    """
    if dms_df is None:
        try:
            dms_df = get_dms_data()
        except FileNotFoundError:
            return pd.DataFrame()

    if dms_df.empty or 'pos' not in dms_df.columns:
        return pd.DataFrame()

    # Group by position
    pos_stats = dms_df.groupby('pos').agg({
        'score': ['mean', 'std', 'min', 'max', 'count'],
        'mutation': 'nunique'
    }).round(3)

    pos_stats.columns = ['mean_score', 'std_score', 'min_score', 'max_score',
                         'n_measurements', 'n_unique_mutations']
    pos_stats = pos_stats.reset_index()

    # Add domain annotation
    pos_stats['in_dbd'] = pos_stats['pos'].apply(
        lambda p: P53_DBD_RANGE[0] <= p <= P53_DBD_RANGE[1]
    )
    pos_stats['is_critical'] = pos_stats['pos'].apply(
        lambda p: p in P53_CRITICAL_POSITIONS
    )

    return pos_stats.sort_values('pos')
