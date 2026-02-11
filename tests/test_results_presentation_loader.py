import unittest

import numpy as np
import pandas as pd

from p53cad.results.presentation import hydrate_display_candidates


class ResultsPresentationLoaderTests(unittest.TestCase):
    def test_hydrate_display_candidates_from_tables(self):
        candidates_df = pd.DataFrame(
            [
                {
                    "candidate_uid": "c1",
                    "profile": "Balanced",
                    "sequence": "AAAA",
                    "score": 1.23,
                    "score_raw": 1.40,
                    "score_calibrated": 1.20,
                    "stability": -0.1,
                    "binding": 7.5,
                    "identity": 96.0,
                    "n_mutations": 2,
                    "mutations_json": '["A1V","A2G"]',
                    "mut_positions_json": "[1,2]",
                    "meets_constraints": True,
                }
            ]
        )
        trajectories_df = pd.DataFrame(
            [
                {
                    "candidate_uid": "c1",
                    "step": 5,
                    "score": 0.8,
                    "mutations_json": '["A1V"]',
                    "mut_positions_json": "[1]",
                }
            ]
        )
        top_df = pd.DataFrame([{"candidate_uid": "c1", "presentation_rank": 1}])

        rows = hydrate_display_candidates(candidates_df, trajectories_df, top_df=top_df)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["candidate_uid"], "c1")
        self.assertEqual(rows[0]["mut_positions"], [1, 2])
        self.assertEqual(len(rows[0]["trajectory"]), 1)

    def test_hydrate_backfills_selection_metadata_when_missing(self):
        candidates_df = pd.DataFrame(
            [
                {
                    "candidate_uid": "c1",
                    "profile": "Balanced",
                    "sequence": "AAAA",
                    "score": 1.20,
                    "mutations_json": '["R175H","N239Y"]',
                    "pass_name": "deep",
                },
                {
                    "candidate_uid": "c2",
                    "profile": "Experimental",
                    "sequence": "BBBB",
                    "score": 1.10,
                    "mutations_json": '["R175H","T284R"]',
                    "pass_name": "deep",
                },
            ]
        )
        rows = hydrate_display_candidates(candidates_df, trajectories_df=pd.DataFrame(), top_df=None)
        self.assertEqual(len(rows), 2)
        for row in rows:
            self.assertTrue(np.isfinite(row["selection_score"]))
            self.assertTrue(np.isfinite(row["diversity_novelty"]))
            self.assertTrue(np.isfinite(row["max_mutation_overlap"]))
            self.assertTrue(str(row["selection_reason"]))


if __name__ == "__main__":
    unittest.main()
