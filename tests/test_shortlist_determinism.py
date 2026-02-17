import unittest

import pandas as pd

from p53cad.results.schema import select_presentation_shortlist


class ShortlistDeterminismTests(unittest.TestCase):
    def test_shortlist_is_deterministic_for_fixed_input(self):
        data = pd.DataFrame(
            [
                {
                    "candidate_uid": f"cand_{i}",
                    "scenario_id": "R175H__gene_therapy" if i < 3 else "R248Q__mrna_therapy",
                    "target_label": "R175H" if i < 3 else "R248Q",
                    "targets_json": "[\"R175H\"]" if i < 3 else "[\"R248Q\"]",
                    "delivery_method": "gene_therapy" if i % 2 == 0 else "mrna_therapy",
                    "profile": "Balanced" if i % 2 == 0 else "Stability-First",
                    "score": 1.0 - i * 0.01,
                    "clinical_score": 80 - i,
                    "uncertainty": 0.02 + i * 0.01,
                    "ood_distance": 0.2 + i * 0.01,
                    "n_mutations": 2,
                    "mutations_json": "[\"R175H\",\"N239Y\"]" if i < 3 else "[\"R248Q\",\"S241H\"]",
                    "sequence": f"SEQ_{i}",
                }
                for i in range(6)
            ]
        )
        s1 = select_presentation_shortlist(data, top_n=4, max_per_target=3)
        s2 = select_presentation_shortlist(data, top_n=4, max_per_target=3)
        self.assertListEqual(s1["candidate_uid"].tolist(), s2["candidate_uid"].tolist())
        self.assertTrue({"selection_score", "diversity_novelty", "max_mutation_overlap", "selection_reason"}.issubset(set(s1.columns)))


if __name__ == "__main__":
    unittest.main()
