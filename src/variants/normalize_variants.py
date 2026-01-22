from __future__ import annotations

import re
from typing import Iterable

import pandas as pd

AA3_TO_AA1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    "SEC": "U",
    "PYL": "O",
    "TER": "*",
}

HGVS_PROTEIN_RE = re.compile(r"p\.[A-Za-z]{1,3}\d+[A-Za-z*]{1,3}(?=$|\s|\)|;|,)")
HGVS_TRANSCRIPT_RE = re.compile(r"([A-Z]{2,4}_\d+\.\d+)")


def _normalize_aa(code: str) -> str:
    if len(code) == 1:
        return code.upper()
    return AA3_TO_AA1.get(code.upper(), "?")


def normalize_protein_change(p_hgvs: str) -> tuple[int, str, str]:
    """Normalize HGVS protein notation into (position, ref, alt)."""
    if not p_hgvs:
        raise ValueError("Missing protein HGVS")

    match = HGVS_PROTEIN_RE.search(p_hgvs)
    if not match:
        raise ValueError(f"Unrecognized protein HGVS: {p_hgvs}")

    token = match.group(0).replace("p.", "")
    token = token.replace("(", "").replace(")", "").strip()

    match_3 = re.match(r"^([A-Za-z]{3})(\d+)([A-Za-z]{3})$", token)
    match_1 = re.match(r"^([A-Za-z\*])(\d+)([A-Za-z\*])$", token)

    if match_3:
        ref3, pos_text, alt3 = match_3.groups()
        ref = _normalize_aa(ref3)
        alt = _normalize_aa(alt3)
    elif match_1:
        ref, pos_text, alt = match_1.groups()
        ref = _normalize_aa(ref)
        alt = _normalize_aa(alt)
    else:
        raise ValueError(f"Unsupported protein HGVS: {p_hgvs}")

    if ref in {"*", "?"} or alt in {"*", "?"}:
        raise ValueError(f"Non-missense protein HGVS: {p_hgvs}")

    pos = int(pos_text)
    return pos, ref, alt


def extract_transcript(hgvs_values: Iterable[str | None]) -> str | None:
    for value in hgvs_values:
        if not value:
            continue
        match = HGVS_TRANSCRIPT_RE.search(value)
        if match:
            return match.group(1)
    return None


def filter_to_canonical_transcript(df: pd.DataFrame, canonical_transcript: str | None = None) -> pd.DataFrame:
    if not canonical_transcript or "transcript" not in df.columns:
        return df.copy()

    canonical_base = canonical_transcript.split(".")[0]
    transcript_series = df["transcript"].fillna("")
    mask = transcript_series.str.startswith(canonical_transcript) | transcript_series.str.startswith(canonical_base)
    return df.loc[mask].copy()


def drop_conflicted_labels(df: pd.DataFrame) -> pd.DataFrame:
    if "clinical_significance" not in df.columns:
        return df.copy()

    mask = ~df["clinical_significance"].fillna("").str.contains("conflict", case=False)
    return df.loc[mask].copy()


def _review_status_score(status: str | None) -> int:
    if not status:
        return 0
    value = status.lower()
    if "practice guideline" in value:
        return 4
    if "reviewed by expert panel" in value:
        return 3
    if "multiple submitters" in value and "no conflicts" in value:
        return 2
    if "criteria provided" in value:
        return 1
    return 0


def dedupe_variants(df: pd.DataFrame) -> pd.DataFrame:
    if not {"pos", "ref", "alt"}.issubset(df.columns):
        return df.copy()

    df = df.copy()
    df["review_score"] = df["review_status"].apply(_review_status_score)
    if "last_evaluated" in df.columns:
        df["last_evaluated_dt"] = pd.to_datetime(df["last_evaluated"], errors="coerce")
    else:
        df["last_evaluated_dt"] = pd.NaT

    df = df.sort_values(["review_score", "last_evaluated_dt"], ascending=[False, False])
    df = df.drop_duplicates(subset=["pos", "ref", "alt"], keep="first")
    return df.drop(columns=["review_score", "last_evaluated_dt"])
