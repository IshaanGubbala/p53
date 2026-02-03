from __future__ import annotations

from typing import Callable, Iterable

MutationSet = tuple[str, ...]


def beam_step(
    current_sets: Iterable[MutationSet],
    generator: Callable[[MutationSet], Iterable[MutationSet]],
    scorer: Callable[[MutationSet], float],
    beam_width: int,
) -> list[tuple[MutationSet, float]]:
    candidates: dict[MutationSet, float] = {}
    for mset in current_sets:
        for candidate in generator(mset):
            if candidate in candidates:
                continue
            candidates[candidate] = scorer(candidate)

    ranked = sorted(candidates.items(), key=lambda item: item[1])
    return ranked[:beam_width]


def run_beam_search(
    seed_mutant: MutationSet,
    generator: Callable[[MutationSet], Iterable[MutationSet]],
    scorer: Callable[[MutationSet], float],
    beam_width: int,
    depth: int,
) -> list[tuple[MutationSet, float]]:
    current = [seed_mutant]
    scored: list[tuple[MutationSet, float]] = [(seed_mutant, scorer(seed_mutant))]

    for _ in range(depth):
        ranked = beam_step(current, generator, scorer, beam_width)
        current = [item[0] for item in ranked]
        scored.extend(ranked)

    seen = set()
    unique: list[tuple[MutationSet, float]] = []
    for mset, score in sorted(scored, key=lambda item: item[1]):
        if mset in seen:
            continue
        seen.add(mset)
        unique.append((mset, score))
    return unique
