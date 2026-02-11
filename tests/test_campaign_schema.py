import unittest

from p53cad.results.schema import (
    BIG8_HOTSPOTS,
    DEFAULT_DELIVERY_METHODS,
    build_scenario_matrix,
    scenarios_to_frame,
)


class CampaignSchemaTests(unittest.TestCase):
    def test_big8_singles_and_pairs_count(self):
        scenarios = build_scenario_matrix(
            hotspots=BIG8_HOTSPOTS,
            delivery_methods=DEFAULT_DELIVERY_METHODS,
            include_pairs=True,
        )
        # 8 singles + C(8,2)=28 pairs => 36 target groups * 3 delivery modes = 108
        self.assertEqual(len(scenarios), 108)

    def test_frame_contains_expected_columns(self):
        scenarios = build_scenario_matrix(include_pairs=False)
        df = scenarios_to_frame(scenarios)
        for col in ["scenario_id", "target_label", "targets", "delivery_method", "n_targets"]:
            self.assertIn(col, df.columns)
        self.assertTrue((df["n_targets"] == 1).all())


if __name__ == "__main__":
    unittest.main()
