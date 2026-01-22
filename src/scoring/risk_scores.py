from __future__ import annotations

import re
from typing import Iterable

MutationSet = tuple[str, ...]


MUTATION_RE = re.compile(r"^([A-Za-z\*])(\d+)([A-Za-z\*])$")


def _parse_mutation(mut: str) -> tuple[int, str, str]:
    # Handle both string and tuple formats
    if isinstance(mut, tuple):
        return mut  # Already in correct format
    
    match = MUTATION_RE.match(str(mut).strip())
    if not match:
        raise ValueError(f"Invalid mutation string: {mut}")
    ref, pos_text, alt = match.groups()
    return int(pos_text), ref.upper(), alt.upper()


def _min_distance(pos: int, protected: Iterable[int], dist_map: dict[tuple[int, int], float]) -> float:
    distances: list[float] = []
    for prot in protected:
        if pos == prot:
            return 0.0
        dist = dist_map.get((pos, prot))
        if dist is None:
            dist = dist_map.get((prot, pos))
        if dist is not None:
            distances.append(dist)
    return min(distances) if distances else float("inf")


def _distance_risk(distance: float, cutoff: float) -> float:
    if distance >= cutoff:
        return 0.0
    if distance <= 0:
        return 1.0
    return max(0.0, 1.0 - distance / cutoff)


def functional_risk(
    mset: MutationSet,
    dist_map: dict[tuple[int, int], float],
    protected: set[int],
    cutoff: float = 8.0,
) -> float:
    if not mset:
        return 0.0
    risks: list[float] = []
    for mut in mset:
        pos, _, _ = _parse_mutation(mut)
        d = _min_distance(pos, protected, dist_map)
        risks.append(_distance_risk(d, cutoff))
    return sum(risks) / len(risks)


def conservation_risk(mset: MutationSet, cons_map: dict[int, float]) -> float:
    if not mset:
        return 0.0
    scores: list[float] = []
    for mut in mset:
        pos, _, _ = _parse_mutation(mut)
        score = cons_map.get(pos)
        if score is None:
            continue
        scores.append(score)
    return sum(scores) / len(scores) if scores else 0.0


def burial_risk(mset: MutationSet, burial_map: dict[int, str]) -> float:
    if not mset:
        return 0.0
    weights = {
        "buried": 0.0,
        "partial": 0.5,
        "exposed": 1.0,
    }
    scores: list[float] = []
    for mut in mset:
        pos, _, _ = _parse_mutation(mut)
        label = burial_map.get(pos)
        if label is None:
            continue
        scores.append(weights.get(label, 0.0))
    return sum(scores) / len(scores) if scores else 0.0


def msa_conservation_risk(mset: MutationSet, msa_cons_map: dict[int, float]) -> float:
    """Calculate risk based on MSA-derived conservation scores.

    Args:
        mset: Mutation set
        msa_cons_map: Position → conservation score (0-1, higher = more conserved)

    Returns:
        Average MSA conservation score (0-1, higher = higher risk)
    """
    if not mset:
        return 0.0

    scores: list[float] = []
    for mut in mset:
        pos, _, _ = _parse_mutation(mut)
        score = msa_cons_map.get(pos)
        if score is None:
            continue
        scores.append(score)

    return sum(scores) / len(scores) if scores else 0.0


def aggregate_risk(weights: dict[str, float], components: dict[str, float]) -> float:
    if not components:
        return 0.0
    total = 0.0
    weight_sum = 0.0
    for name, value in components.items():
        weight = weights.get(name, 1.0)
        total += weight * value
        weight_sum += weight
    return total / weight_sum if weight_sum else 0.0
