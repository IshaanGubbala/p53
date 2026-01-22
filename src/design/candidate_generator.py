from __future__ import annotations

import re
from typing import Iterable

from src.design.mutation_sets import format_mutation

MutationSet = tuple[str, ...]


CONSERVATIVE_MAP = {
    "A": ["V", "L", "S", "T", "G"],
    "R": ["K", "H", "Q"],
    "N": ["D", "Q", "S"],
    "D": ["E", "N"],
    "C": ["S", "A"],
    "Q": ["E", "N", "K"],
    "E": ["D", "Q"],
    "G": ["A", "S"],
    "H": ["R", "K", "N"],
    "I": ["V", "L", "M"],
    "L": ["I", "V", "M"],
    "K": ["R", "H", "Q"],
    "M": ["L", "I", "V"],
    "F": ["Y", "W", "L"],
    "P": ["A", "S"],
    "S": ["T", "A", "N"],
    "T": ["S", "A", "V"],
    "W": ["F", "Y"],
    "Y": ["F", "W"],
    "V": ["I", "L", "M"],
}


_MUT_RE = re.compile(r"^([A-Za-z\*])(\d+)([A-Za-z\*])$")


def select_design_positions(
    burial_map: dict[int, str],
    protected: set[int],
    cons_map: dict[int, float] | None = None,
    top_n: int = 80,
    max_cons: float = 0.8,
    allow_exposed: bool = False,
) -> list[int]:
    cons_map = cons_map or {}

    scores: list[tuple[float, int]] = []
    for pos, label in burial_map.items():
        if pos in protected:
            continue
        if label == "exposed" and not allow_exposed:
            continue
        cons_score = cons_map.get(pos, 0.0)
        if cons_score > max_cons:
            continue
        burial_score = {"buried": 2.0, "partial": 1.0, "exposed": 0.0}.get(label, 0.0)
        score = burial_score - cons_score
        scores.append((score, pos))

    scores.sort(reverse=True)
    return [pos for _, pos in scores[:top_n]]


def propose_substitutions(pos: int, ref_aa: str) -> list[str]:
    ref = ref_aa.upper()
    if ref not in CONSERVATIVE_MAP:
        return []
    return CONSERVATIVE_MAP[ref]


def generate_single_mutants(positions: Iterable[int], sequence: str) -> list[MutationSet]:
    muts: list[MutationSet] = []
    for pos in positions:
        if pos < 1 or pos > len(sequence):
            continue
        ref = sequence[pos - 1]
        for alt in propose_substitutions(pos, ref):
            if alt == ref:
                continue
            muts.append((format_mutation(pos, ref, alt),))
    return muts


def _mutation_pos(mut: str) -> int:
    match = _MUT_RE.match(mut.strip())
    if not match:
        raise ValueError(f"Invalid mutation: {mut}")
    return int(match.group(2))


def expand_to_pairs(best_singles: list[MutationSet], max_pairs: int) -> list[MutationSet]:
    pairs: list[MutationSet] = []
    for i in range(len(best_singles)):
        pos_i = _mutation_pos(best_singles[i][0])
        for j in range(i + 1, len(best_singles)):
            pos_j = _mutation_pos(best_singles[j][0])
            if pos_i == pos_j:
                continue
            pair = tuple(sorted({best_singles[i][0], best_singles[j][0]}))
            pairs.append(pair)
            if len(pairs) >= max_pairs:
                return pairs
    return pairs


def expand_to_triples(best_pairs: list[MutationSet], max_triples: int) -> list[MutationSet]:
    triples: list[MutationSet] = []
    for i in range(len(best_pairs)):
        muts_i = set(best_pairs[i])
        positions_i = {_mutation_pos(mut) for mut in muts_i}
        for j in range(i + 1, len(best_pairs)):
            muts_j = set(best_pairs[j])
            positions_j = {_mutation_pos(mut) for mut in muts_j}
            if positions_i & positions_j:
                continue
            combined = tuple(sorted(muts_i | muts_j))
            if len(combined) != 3:
                continue
            triples.append(combined)
            if len(triples) >= max_triples:
                return triples
    return triples
