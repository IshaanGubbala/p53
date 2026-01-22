from __future__ import annotations

from pathlib import Path

import pandas as pd


def load_uniprot_sequence(fasta_path: Path) -> str:
    """Load a UniProt FASTA sequence (single record)."""
    sequence: list[str] = []
    with fasta_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                continue
            sequence.append(line.strip())
    return "".join(sequence)


def validate_variant_ref(seq: str, pos: int, ref_aa: str) -> bool:
    """Validate reference amino acid at 1-based position."""
    if pos < 1 or pos > len(seq):
        return False
    return seq[pos - 1].upper() == ref_aa.upper()


def map_and_filter(df: pd.DataFrame, seq: str) -> pd.DataFrame:
    """Drop variants with mismatched reference residues."""
    if not {"pos", "ref"}.issubset(df.columns):
        return df.copy()
    if df.empty:
        return df.copy()

    pos_series = pd.to_numeric(df["pos"], errors="coerce")
    ref_series = df["ref"].fillna("").astype(str).str.upper()

    seq_len = len(seq)
    ref_lookup: list[str] = []
    for pos in pos_series:
        if pd.isna(pos):
            ref_lookup.append("")
            continue
        pos_int = int(pos)
        if pos_int < 1 or pos_int > seq_len:
            ref_lookup.append("")
            continue
        ref_lookup.append(seq[pos_int - 1].upper())

    mask = ref_series == pd.Series(ref_lookup, index=df.index)
    return df.loc[mask].copy()


def restrict_to_domain(df: pd.DataFrame, start: int, end: int) -> pd.DataFrame:
    """Restrict variants to a 1-based inclusive domain range."""
    if "pos" not in df.columns:
        return df.copy()

    mask = (df["pos"] >= start) & (df["pos"] <= end)
    return df.loc[mask].copy()
