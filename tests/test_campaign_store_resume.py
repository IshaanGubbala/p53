import tempfile
import unittest
from pathlib import Path

import pandas as pd

from p53cad.results.store import CampaignStore


class CampaignStoreResumeTests(unittest.TestCase):
    def test_checkpoint_resume_roundtrip(self):
        with tempfile.TemporaryDirectory() as tmp:
            store = CampaignStore(base_dir=Path(tmp) / "campaigns")
            run_id = "campaign_test"
            store.init_run(run_id, config={"seed": 42}, runtime_caps={}, resume=False)
            store.append_scenario_checkpoint(
                run_id,
                {
                    "scenario_id": "R175H__gene_therapy",
                    "pass_name": "screen",
                    "status": "done",
                },
            )
            ckpt = store.load_scenario_checkpoints(run_id)
            self.assertFalse(ckpt.empty)
            self.assertEqual(ckpt.iloc[0]["status"], "done")

    def test_write_and_load_bundle_tables(self):
        with tempfile.TemporaryDirectory() as tmp:
            store = CampaignStore(base_dir=Path(tmp) / "campaigns")
            run_id = "campaign_test"
            store.init_run(run_id, config={"seed": 42}, runtime_caps={}, resume=False)
            store.write_tables(
                run_id,
                scenarios_df=pd.DataFrame([{"scenario_id": "a"}]),
                candidates_df=pd.DataFrame([{"candidate_uid": "c1", "score": 1.0}]),
                trajectories_df=pd.DataFrame([{"candidate_uid": "c1", "step": 5}]),
                clinical_df=pd.DataFrame([{"candidate_uid": "c1", "clinical_score": 75.0}]),
            )
            store.write_top30(run_id, pd.DataFrame([{"candidate_uid": "c1", "presentation_rank": 1}]))
            bundle = store.load_run_bundle(run_id)
            self.assertIn("candidates", bundle)
            self.assertFalse(bundle["candidates"].empty)
            self.assertFalse(bundle["top30"].empty)


if __name__ == "__main__":
    unittest.main()
