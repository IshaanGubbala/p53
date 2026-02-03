from __future__ import annotations

from typing import Optional


def consensus_ddg(ddg_foldx: float, ddg_other: Optional[float] = None) -> float:
    if ddg_other is None:
        return ddg_foldx
    return (ddg_foldx + ddg_other) / 2.0


def rank_stability(ddg: float) -> float:
    """Lower ddG indicates higher stability; return a score for ranking."""
    return -ddg
