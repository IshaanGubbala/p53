"""Tests for NSGA-II Pareto optimization module.

Covers dominance relation, non-dominated sorting, crowding distance,
front queries, and batch operations.
"""

import math
import unittest

from p53cad.engine.pareto import ParetoFront, ParetoSolution


class TestDominance(unittest.TestCase):
    """Test the Pareto dominance relation."""

    def _front(self):
        return ParetoFront(objectives=["a", "b"])

    def test_dominates_strictly_better(self):
        pf = self._front()
        a = ParetoSolution("A", {"a": 2.0, "b": 2.0})
        b = ParetoSolution("B", {"a": 1.0, "b": 1.0})
        self.assertTrue(pf._dominates(a, b))
        self.assertFalse(pf._dominates(b, a))

    def test_no_dominance_when_equal(self):
        pf = self._front()
        a = ParetoSolution("A", {"a": 1.0, "b": 1.0})
        b = ParetoSolution("B", {"a": 1.0, "b": 1.0})
        self.assertFalse(pf._dominates(a, b))
        self.assertFalse(pf._dominates(b, a))

    def test_no_dominance_tradeoff(self):
        pf = self._front()
        a = ParetoSolution("A", {"a": 2.0, "b": 0.5})
        b = ParetoSolution("B", {"a": 0.5, "b": 2.0})
        self.assertFalse(pf._dominates(a, b))
        self.assertFalse(pf._dominates(b, a))

    def test_dominates_one_equal_one_better(self):
        pf = self._front()
        a = ParetoSolution("A", {"a": 2.0, "b": 1.0})
        b = ParetoSolution("B", {"a": 1.0, "b": 1.0})
        self.assertTrue(pf._dominates(a, b))


class TestNSGAIIRanking(unittest.TestCase):
    """Test non-dominated sorting assigns correct Pareto ranks."""

    def test_single_solution_rank_1(self):
        pf = ParetoFront(objectives=["x", "y"])
        pf.add(ParetoSolution("only", {"x": 1.0, "y": 1.0}))
        self.assertEqual(pf.solutions[0].pareto_rank, 1)

    def test_two_nondominated_both_rank_1(self):
        pf = ParetoFront(objectives=["x", "y"])
        pf.add(ParetoSolution("A", {"x": 2.0, "y": 0.5}))
        pf.add(ParetoSolution("B", {"x": 0.5, "y": 2.0}))
        sols = pf.solutions
        self.assertEqual(sols[0].pareto_rank, 1)
        self.assertEqual(sols[1].pareto_rank, 1)

    def test_dominated_gets_rank_2(self):
        pf = ParetoFront(objectives=["x", "y"])
        pf.add(ParetoSolution("dominant", {"x": 3.0, "y": 3.0}))
        pf.add(ParetoSolution("dominated", {"x": 1.0, "y": 1.0}))
        sols = {s.candidate_uid: s for s in pf.solutions}
        self.assertEqual(sols["dominant"].pareto_rank, 1)
        self.assertEqual(sols["dominated"].pareto_rank, 2)

    def test_three_rank_layers(self):
        pf = ParetoFront(objectives=["x", "y"])
        pf.add(ParetoSolution("top", {"x": 5.0, "y": 5.0}))
        pf.add(ParetoSolution("mid", {"x": 3.0, "y": 3.0}))
        pf.add(ParetoSolution("bot", {"x": 1.0, "y": 1.0}))
        sols = {s.candidate_uid: s for s in pf.solutions}
        self.assertEqual(sols["top"].pareto_rank, 1)
        self.assertEqual(sols["mid"].pareto_rank, 2)
        self.assertEqual(sols["bot"].pareto_rank, 3)

    def test_batch_add(self):
        pf = ParetoFront(objectives=["x"])
        pf.add_batch([
            ParetoSolution("a", {"x": 1.0}),
            ParetoSolution("b", {"x": 2.0}),
            ParetoSolution("c", {"x": 3.0}),
        ])
        sols = {s.candidate_uid: s for s in pf.solutions}
        self.assertEqual(sols["c"].pareto_rank, 1)
        self.assertEqual(sols["b"].pareto_rank, 2)
        self.assertEqual(sols["a"].pareto_rank, 3)


class TestCrowdingDistance(unittest.TestCase):
    """Test crowding distance computation."""

    def test_boundary_solutions_infinite_distance(self):
        pf = ParetoFront(objectives=["x", "y"])
        # 4 non-dominated solutions on rank 1 front
        pf.add(ParetoSolution("A", {"x": 0.0, "y": 5.0}))
        pf.add(ParetoSolution("B", {"x": 2.0, "y": 3.0}))
        pf.add(ParetoSolution("C", {"x": 3.0, "y": 2.0}))
        pf.add(ParetoSolution("D", {"x": 5.0, "y": 0.0}))
        sols = {s.candidate_uid: s for s in pf.solutions}
        # At least A and D should be boundary in at least one objective
        self.assertTrue(math.isinf(sols["A"].crowding_distance))
        self.assertTrue(math.isinf(sols["D"].crowding_distance))

    def test_two_solutions_infinite_distance(self):
        pf = ParetoFront(objectives=["x"])
        pf.add(ParetoSolution("A", {"x": 1.0}))
        pf.add(ParetoSolution("B", {"x": 2.0}))
        sols = pf.solutions
        self.assertTrue(math.isinf(sols[0].crowding_distance))
        self.assertTrue(math.isinf(sols[1].crowding_distance))

    def test_interior_solution_finite_distance(self):
        pf = ParetoFront(objectives=["x", "y"])
        pf.add(ParetoSolution("lo", {"x": 0.0, "y": 10.0}))
        pf.add(ParetoSolution("mid", {"x": 5.0, "y": 5.0}))
        pf.add(ParetoSolution("hi", {"x": 10.0, "y": 0.0}))
        sols = {s.candidate_uid: s for s in pf.solutions}
        self.assertTrue(math.isfinite(sols["mid"].crowding_distance))
        self.assertGreater(sols["mid"].crowding_distance, 0)


class TestFrontQueries(unittest.TestCase):
    """Test get_front and best_by_crowding queries."""

    def test_get_front_rank_1(self):
        pf = ParetoFront(objectives=["x", "y"])
        pf.add(ParetoSolution("top1", {"x": 5.0, "y": 1.0}))
        pf.add(ParetoSolution("top2", {"x": 1.0, "y": 5.0}))
        pf.add(ParetoSolution("low", {"x": 0.5, "y": 0.5}))
        rank1 = pf.get_front(rank=1)
        uids = {s.candidate_uid for s in rank1}
        self.assertIn("top1", uids)
        self.assertIn("top2", uids)
        self.assertNotIn("low", uids)

    def test_best_by_crowding_respects_top_n(self):
        pf = ParetoFront(objectives=["x"])
        for i in range(20):
            pf.add(ParetoSolution(f"s{i}", {"x": float(i)}))
        best5 = pf.best_by_crowding(top_n=5)
        self.assertEqual(len(best5), 5)
        # Highest rank should come first (rank 1 = x=19)
        self.assertEqual(best5[0].candidate_uid, "s19")

    def test_summary_statistics(self):
        pf = ParetoFront(objectives=["x", "y"])
        pf.add(ParetoSolution("A", {"x": 3.0, "y": 3.0}))
        pf.add(ParetoSolution("B", {"x": 1.0, "y": 1.0}))
        summary = pf.summary()
        self.assertEqual(summary["n_solutions"], 2)
        self.assertEqual(summary["n_fronts"], 2)
        self.assertEqual(summary["rank_1_count"], 1)

    def test_empty_front(self):
        pf = ParetoFront()
        self.assertEqual(len(pf), 0)
        self.assertEqual(pf.summary(), {"n_solutions": 0, "n_fronts": 0})

    def test_five_objective_ranking(self):
        """Test with the real 5-objective setup used in campaigns."""
        pf = ParetoFront()  # defaults to OBJECTIVES
        pf.add(ParetoSolution("best", {
            "oracle_score": 1.0, "pll": -50.0, "dms_quality": 2.0,
            "identity": 98.0, "stability": 5.0,
        }))
        pf.add(ParetoSolution("worst", {
            "oracle_score": 0.1, "pll": -200.0, "dms_quality": -1.0,
            "identity": 80.0, "stability": -2.0,
        }))
        pf.add(ParetoSolution("tradeoff", {
            "oracle_score": 0.8, "pll": -30.0, "dms_quality": 1.0,
            "identity": 85.0, "stability": 8.0,
        }))
        sols = {s.candidate_uid: s for s in pf.solutions}
        # best and tradeoff are both rank 1 (neither dominates the other:
        # best has higher oracle/dms/identity, tradeoff has higher pll/stability)
        self.assertEqual(sols["best"].pareto_rank, 1)
        self.assertEqual(sols["tradeoff"].pareto_rank, 1)
        # worst is dominated by best (and tradeoff), so rank 2
        self.assertEqual(sols["worst"].pareto_rank, 2)


if __name__ == "__main__":
    unittest.main()
