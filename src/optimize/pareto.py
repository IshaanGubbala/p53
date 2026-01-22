from __future__ import annotations

from typing import Iterable


def dominates(a: dict, b: dict, objectives: Iterable[str]) -> bool:
    better_or_equal = True
    strictly_better = False
    for key in objectives:
        aval = a.get(key)
        bval = b.get(key)
        if aval is None or bval is None:
            return False
        if aval > bval:
            better_or_equal = False
            break
        if aval < bval:
            strictly_better = True
    return better_or_equal and strictly_better


def pareto_front(items: list[dict], objectives: Iterable[str]) -> list[dict]:
    front: list[dict] = []
    for item in items:
        dominated = False
        for other in items:
            if other is item:
                continue
            if dominates(other, item, objectives):
                dominated = True
                break
        if not dominated:
            front.append(item)
    return front


def crowding_distance(front: list[dict], objectives: Iterable[str]) -> list[float]:
    if not front:
        return []
    distances = [0.0 for _ in front]
    objective_list = list(objectives)

    for key in objective_list:
        sorted_idx = sorted(range(len(front)), key=lambda i: front[i].get(key, 0.0))
        distances[sorted_idx[0]] = float("inf")
        distances[sorted_idx[-1]] = float("inf")
        values = [front[i].get(key, 0.0) for i in sorted_idx]
        min_val = min(values)
        max_val = max(values)
        if max_val == min_val:
            continue
        for i in range(1, len(sorted_idx) - 1):
            prev_val = values[i - 1]
            next_val = values[i + 1]
            distances[sorted_idx[i]] += (next_val - prev_val) / (max_val - min_val)

    return distances
