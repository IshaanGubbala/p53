"""Tests for autoregressive (Gibbs-like) mutation sampling.

Validates that _run_autoregressive_trial produces valid candidate dicts
with expected fields and constraints.

NOTE: Requires PyTorch, ESM-2 model, and full p53cad pipeline.
"""

import os
import shutil
import tempfile
import unittest

os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

try:
    import torch
    from p53cad.core.runtime import bootstrap_runtime

    bootstrap_runtime()
    from p53cad.data.dms import P53_WT
    from p53cad.engine.campaign import CampaignRunner, ScenarioRuntime
    from p53cad.results.store import CampaignStore

    HAS_DEPS = True
except ImportError:
    HAS_DEPS = False


@unittest.skipUnless(HAS_DEPS, "Requires torch, transformers, p53cad")
class TestAutoRegressiveTrial(unittest.TestCase):
    """Test the autoregressive mutation sampling method."""

    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.mkdtemp(prefix="p53cad_ar_test_")
        cls.store = CampaignStore(base_dir=cls.tmpdir)
        cls.runner = CampaignRunner(store=cls.store)

        # Ensure embedder + oracle are loaded
        cls.runner._load_models()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmpdir, ignore_errors=True)

    def test_produces_candidate_dict(self):
        """AR trial returns a dict with required candidate fields."""
        from p53cad.data.dms import apply_mutation

        target_seq = apply_mutation(P53_WT, "R175H")
        locked_indices = [174]  # 0-indexed position of R175H

        scenario = ScenarioRuntime(
            scenario_id="R175H__gene_therapy",
            target_label="R175H",
            targets=["R175H"],
            delivery_method="gene_therapy",
        )

        cand, traj = self.runner._run_autoregressive_trial(
            scenario=scenario,
            target_seq=target_seq,
            locked_indices=locked_indices,
            max_mutations=5,
            pass_name="ar_test",
            trial_idx=0,
            trial_seed=42,
            baseline_score=0.5,
        )

        # Required fields
        self.assertIn("candidate_uid", cand)
        self.assertIn("sequence", cand)
        self.assertIn("score", cand)
        self.assertIn("mutations_json", cand)
        self.assertIn("n_mutations", cand)
        self.assertIn("identity", cand)

    def test_sequence_length_preserved(self):
        """AR trial produces sequence of correct length (393 for p53)."""
        from p53cad.data.dms import apply_mutation

        target_seq = apply_mutation(P53_WT, "R175H")
        scenario = ScenarioRuntime(
            scenario_id="R175H__gene_therapy",
            target_label="R175H",
            targets=["R175H"],
            delivery_method="gene_therapy",
        )

        cand, traj = self.runner._run_autoregressive_trial(
            scenario=scenario,
            target_seq=target_seq,
            locked_indices=[174],
            max_mutations=3,
            pass_name="ar_len",
            trial_idx=0,
            trial_seed=99,
            baseline_score=0.5,
        )

        self.assertEqual(len(cand["sequence"]), len(P53_WT))

    def test_locked_positions_preserved(self):
        """AR trial preserves the cancer mutation at locked positions."""
        from p53cad.data.dms import apply_mutation

        target_seq = apply_mutation(P53_WT, "R175H")
        scenario = ScenarioRuntime(
            scenario_id="R175H__gene_therapy",
            target_label="R175H",
            targets=["R175H"],
            delivery_method="gene_therapy",
        )

        cand, traj = self.runner._run_autoregressive_trial(
            scenario=scenario,
            target_seq=target_seq,
            locked_indices=[174],
            max_mutations=5,
            pass_name="ar_lock",
            trial_idx=0,
            trial_seed=123,
            baseline_score=0.5,
        )

        # Position 174 (0-indexed) = residue 175, should be H (histidine)
        self.assertEqual(cand["sequence"][174], "H")

    def test_score_is_finite(self):
        """AR trial candidate score should be a finite number."""
        import math
        from p53cad.data.dms import apply_mutation

        target_seq = apply_mutation(P53_WT, "R273H")
        scenario = ScenarioRuntime(
            scenario_id="R273H__mrna_therapy",
            target_label="R273H",
            targets=["R273H"],
            delivery_method="mrna_therapy",
        )

        cand, traj = self.runner._run_autoregressive_trial(
            scenario=scenario,
            target_seq=target_seq,
            locked_indices=[272],
            max_mutations=4,
            pass_name="ar_finite",
            trial_idx=0,
            trial_seed=7,
            baseline_score=0.5,
        )

        self.assertFalse(math.isnan(cand["score"]))
        self.assertFalse(math.isinf(cand["score"]))


if __name__ == "__main__":
    unittest.main()
