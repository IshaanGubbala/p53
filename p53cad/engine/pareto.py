"""Multi-objective Pareto optimization for p53 rescue candidate ranking.

Implements NSGA-II non-dominated sorting and crowding distance for ranking
candidates across multiple objectives (oracle score, PLL, DMS quality,
identity, stability).  All objectives are expressed as values to MAXIMIZE.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence

import numpy as np


@dataclass
class ParetoSolution:
    """A single candidate in multi-objective space."""

    candidate_uid: str
    objectives: Dict[str, float]  # all values to MAXIMIZE
    pareto_rank: int = 0
    crowding_distance: float = 0.0
    metadata: Dict[str, Any] = field(default_factory=dict)


class ParetoFront:
    """NSGA-II non-dominated sorting and crowding distance computation.

    Usage::

        front = ParetoFront()
        for cand in candidates:
            front.add(ParetoSolution(
                candidate_uid=cand["candidate_uid"],
                objectives={
                    "oracle_score": cand["score"],
                    "pll": cand.get("pll", 0.0),
                    "dms_quality": -cand.get("rescue_dms_mean", 0.0),  # negate: lower Z = better
                    "identity": cand.get("identity", 100.0),
                    "stability": cand.get("stability", 0.0),
                },
            ))
        rank1 = front.get_front(rank=1)
        best30 = front.best_by_crowding(top_n=30)
    """

    OBJECTIVES = ["oracle_score", "pll", "dms_quality", "identity", "stability"]

    def __init__(self, objectives: Optional[Sequence[str]] = None):
        self._objectives = list(objectives or self.OBJECTIVES)
        self._solutions: List[ParetoSolution] = []
        self._dirty = True

    def add(self, solution: ParetoSolution) -> None:
        """Add a solution and mark ranks as stale."""
        self._solutions.append(solution)
        self._dirty = True

    def add_batch(self, solutions: Sequence[ParetoSolution]) -> None:
        """Add multiple solutions at once."""
        self._solutions.extend(solutions)
        self._dirty = True

    @property
    def solutions(self) -> List[ParetoSolution]:
        if self._dirty:
            self._recompute_ranks()
        return list(self._solutions)

    def _dominates(self, a: ParetoSolution, b: ParetoSolution) -> bool:
        """Return True if solution *a* Pareto-dominates *b*.

        *a* dominates *b* iff a is >= b on all objectives and strictly > on at least one.
        """
        at_least_one_better = False
        for obj in self._objectives:
            va = a.objectives.get(obj, 0.0)
            vb = b.objectives.get(obj, 0.0)
            if va < vb:
                return False
            if va > vb:
                at_least_one_better = True
        return at_least_one_better

    def _recompute_ranks(self) -> None:
        """Assign NSGA-II Pareto ranks via non-dominated sorting."""
        n = len(self._solutions)
        if n == 0:
            self._dirty = False
            return

        # domination_count[i] = number of solutions that dominate i
        domination_count = [0] * n
        # dominated_set[i] = indices of solutions that i dominates
        dominated_set: List[List[int]] = [[] for _ in range(n)]

        for i in range(n):
            for j in range(i + 1, n):
                if self._dominates(self._solutions[i], self._solutions[j]):
                    dominated_set[i].append(j)
                    domination_count[j] += 1
                elif self._dominates(self._solutions[j], self._solutions[i]):
                    dominated_set[j].append(i)
                    domination_count[i] += 1

        # Build fronts
        current_front: List[int] = []
        for i in range(n):
            if domination_count[i] == 0:
                self._solutions[i].pareto_rank = 1
                current_front.append(i)

        rank = 1
        while current_front:
            next_front: List[int] = []
            for i in current_front:
                for j in dominated_set[i]:
                    domination_count[j] -= 1
                    if domination_count[j] == 0:
                        self._solutions[j].pareto_rank = rank + 1
                        next_front.append(j)
            rank += 1
            current_front = next_front

        self._compute_crowding_distances()
        self._dirty = False

    def _compute_crowding_distances(self) -> None:
        """Compute crowding distance within each Pareto rank."""
        if not self._solutions:
            return

        max_rank = max(s.pareto_rank for s in self._solutions)

        for rank in range(1, max_rank + 1):
            indices = [i for i, s in enumerate(self._solutions) if s.pareto_rank == rank]
            if len(indices) <= 2:
                for i in indices:
                    self._solutions[i].crowding_distance = float("inf")
                continue

            # Initialize distances
            distances = {i: 0.0 for i in indices}

            for obj in self._objectives:
                # Sort indices by this objective
                sorted_idx = sorted(indices, key=lambda i: self._solutions[i].objectives.get(obj, 0.0))

                # Boundary solutions get infinite distance
                distances[sorted_idx[0]] = float("inf")
                distances[sorted_idx[-1]] = float("inf")

                # Objective range for normalization
                obj_min = self._solutions[sorted_idx[0]].objectives.get(obj, 0.0)
                obj_max = self._solutions[sorted_idx[-1]].objectives.get(obj, 0.0)
                obj_range = obj_max - obj_min
                if obj_range < 1e-12:
                    continue

                for k in range(1, len(sorted_idx) - 1):
                    prev_val = self._solutions[sorted_idx[k - 1]].objectives.get(obj, 0.0)
                    next_val = self._solutions[sorted_idx[k + 1]].objectives.get(obj, 0.0)
                    distances[sorted_idx[k]] += (next_val - prev_val) / obj_range

            for i in indices:
                self._solutions[i].crowding_distance = distances[i]

    def get_front(self, rank: int = 1) -> List[ParetoSolution]:
        """Return all solutions with the given Pareto rank."""
        if self._dirty:
            self._recompute_ranks()
        return [s for s in self._solutions if s.pareto_rank == rank]

    def best_by_crowding(self, top_n: int = 30) -> List[ParetoSolution]:
        """Return top-N solutions sorted by (rank asc, crowding_distance desc)."""
        if self._dirty:
            self._recompute_ranks()

        sorted_solutions = sorted(
            self._solutions,
            key=lambda s: (s.pareto_rank, -s.crowding_distance),
        )
        return sorted_solutions[:top_n]

    def __len__(self) -> int:
        return len(self._solutions)

    def summary(self) -> Dict[str, Any]:
        """Return summary statistics of the Pareto front."""
        if self._dirty:
            self._recompute_ranks()
        if not self._solutions:
            return {"n_solutions": 0, "n_fronts": 0}

        ranks = [s.pareto_rank for s in self._solutions]
        return {
            "n_solutions": len(self._solutions),
            "n_fronts": max(ranks),
            "rank_1_count": sum(1 for r in ranks if r == 1),
            "rank_distribution": {r: sum(1 for rr in ranks if rr == r) for r in sorted(set(ranks))},
        }
