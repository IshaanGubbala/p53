import unittest

import pandas as pd

from p53cad.engine.campaign import CampaignRunner


class CampaignSummaryTests(unittest.TestCase):
    def test_summary_reports_unique_scenario_and_pass_breakdown(self):
        runner = CampaignRunner.__new__(CampaignRunner)
        scenario_df = pd.DataFrame(
            [
                {"scenario_id": "R175H__gene_therapy", "pass_name": "screen"},
                {"scenario_id": "R175H__gene_therapy", "pass_name": "deep"},
                {"scenario_id": "R248Q__mrna_therapy", "pass_name": "screen"},
            ]
        )
        candidate_df = pd.DataFrame([{"candidate_uid": "c1", "score": 1.0, "identity": 99.0, "n_mutations": 2}])
        top_df = pd.DataFrame([{"presentation_rank": 1, "target_label": "R175H", "delivery_method": "gene_therapy", "score": 1.0}])

        summary = runner._build_summary("run_x", scenario_df, candidate_df, top_df)
        self.assertIn("- Scenarios (unique): 2", summary)
        self.assertIn("- Scenario pass rows:", summary)


if __name__ == "__main__":
    unittest.main()
