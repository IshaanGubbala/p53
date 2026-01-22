from __future__ import annotations

import re
from typing import Iterable


def _mutation_pos(mut: str) -> int:
    match = re.match(r"^([A-Za-z\*])(\d+)([A-Za-z\*])$", mut.strip())
    if not match:
        raise ValueError(f"Invalid mutation: {mut}")
    return int(match.group(2))


def filter_by_protected_distance(
    mset: Iterable[str],
    dist_map: dict[tuple[int, int], float],
    protected: set[int],
    min_angstrom: float,
) -> bool:
    for mut in mset:
        pos = _mutation_pos(mut)
        if pos in protected:
            return False
        min_dist = float("inf")
        for prot in protected:
            dist = dist_map.get((pos, prot))
            if dist is None:
                dist = dist_map.get((prot, pos))
            if dist is not None:
                min_dist = min(min_dist, dist)
        if min_dist < min_angstrom:
            return False
    return True


def filter_by_conservation(
    mset: Iterable[str],
    cons_map: dict[int, float],
    max_cons: float,
) -> bool:
    if not cons_map:
        return True
    for mut in mset:
        pos = _mutation_pos(mut)
        score = cons_map.get(pos)
        if score is None:
            continue
        if score > max_cons:
            return False
    return True


def filter_by_surface_patch(
    mset: Iterable[str],
    burial_map: dict[int, str],
    allow_exposed: bool = False,
) -> bool:
    for mut in mset:
        pos = _mutation_pos(mut)
        label = burial_map.get(pos)
        if label is None:
            continue
        if label == "exposed" and not allow_exposed:
            return False
    return True


def apply_all_filters(
    candidates: Iterable[tuple[str, ...]],
    dist_map: dict[tuple[int, int], float],
    protected: set[int],
    burial_map: dict[int, str],
    min_angstrom: float,
    cons_map: dict[int, float] | None = None,
    max_cons: float = 0.8,
    allow_exposed: bool = False,
) -> list[tuple[str, ...]]:
    cons_map = cons_map or {}
    filtered: list[tuple[str, ...]] = []
    for mset in candidates:
        if not filter_by_protected_distance(mset, dist_map, protected, min_angstrom):
            continue
        if not filter_by_conservation(mset, cons_map, max_cons):
            continue
        if not filter_by_surface_patch(mset, burial_map, allow_exposed=allow_exposed):
            continue
        filtered.append(mset)
    return filtered
