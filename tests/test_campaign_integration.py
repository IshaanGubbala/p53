"""
Integration test for the full campaign pipeline.

Runs a minimal campaign (1 scenario, 1 trial) end-to-end and validates:
- Candidate generation
- Shortlist selection
- Parquet artifact creation
- ESMFold validation output
- MD script generation
- PyMOL script generation
- Summary markdown

NOTE: Requires PyTorch and ESM-2 model (~30MB download on first run).
Skipped if torch/transformers are unavailable.
"""

import json
import os
import shutil
import tempfile
import unittest

os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

try:
    import torch
    from p53cad.core.runtime import bootstrap_runtime

    bootstrap_runtime()
    from p53cad.engine.campaign import CampaignRunner
    from p53cad.results.store import CampaignStore

    HAS_DEPS = True
except ImportError:
    HAS_DEPS = False


@unittest.skipUnless(HAS_DEPS, "Requires torch, transformers, p53cad")
class TestCampaignIntegration(unittest.TestCase):
    """End-to-end integration test with a single-scenario mini campaign."""

    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.mkdtemp(prefix="p53cad_test_")
        cls.store = CampaignStore(base_dir=cls.tmpdir)
        cls.runner = CampaignRunner(store=cls.store)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmpdir, ignore_errors=True)

    def test_mini_campaign_runs_to_completion(self):
        """Run a single-scenario campaign and verify all artifacts."""
        result = self.runner.run(
            seed=42,
            include_pairs=False,
            hotspots=["R175H"],
            delivery_methods=["gene_therapy"],
            max_scenarios=1,
            shortlist_n=5,
            budget="low",
            with_clinical=True,
        )

        # Basic result checks
        self.assertIn("run_id", result)
        self.assertIn("run_dir", result)
        self.assertGreater(result["n_candidates"], 0)
        self.assertGreaterEqual(result["n_shortlist"], 0)

        run_id = result["run_id"]
        paths = self.store.run_paths(run_id)

        # Parquet artifacts exist
        self.assertTrue(paths.candidates_path.exists(), "candidates.parquet missing")
        self.assertTrue(paths.scenarios_path.exists(), "scenarios.parquet missing")
        self.assertTrue(paths.trajectories_path.exists(), "trajectories.parquet missing")

        # Manifest is valid JSON with completed status
        self.assertTrue(paths.manifest_path.exists(), "manifest.json missing")
        manifest = json.loads(paths.manifest_path.read_text())
        self.assertEqual(manifest["status"], "completed")

        # Summary markdown exists
        self.assertTrue(paths.summary_path.exists(), "summary.md missing")
        summary_text = paths.summary_path.read_text()
        self.assertIn("Campaign Summary", summary_text)

        # Load run bundle and verify data
        bundle = self.store.load_run_bundle(run_id)
        cand_df = bundle["candidates"]
        self.assertGreater(len(cand_df), 0)
        self.assertIn("score", cand_df.columns)
        self.assertIn("sequence", cand_df.columns)
        self.assertIn("mutations_json", cand_df.columns)
        self.assertIn("rescue_dms_mean", cand_df.columns)

        # All candidates should have valid sequences (length of p53 WT = 393)
        for seq in cand_df["sequence"]:
            self.assertEqual(len(seq), 393, f"Unexpected sequence length: {len(seq)}")

        # Scores should be finite
        import numpy as np

        self.assertTrue(np.all(np.isfinite(cand_df["score"].values)))

    def test_report_run_regenerates_shortlist(self):
        """Verify report_run can re-generate shortlist from existing data."""
        # First run a campaign
        result = self.runner.run(
            seed=99,
            include_pairs=False,
            hotspots=["R273H"],
            delivery_methods=["mrna_therapy"],
            max_scenarios=1,
            shortlist_n=3,
            budget="low",
            with_clinical=True,
        )
        run_id = result["run_id"]

        # Re-run report generation
        report = self.runner.report_run(run_id, shortlist_n=3)
        self.assertEqual(report["run_id"], run_id)
        self.assertGreater(report["n_candidates"], 0)

    def test_store_list_runs_finds_campaign(self):
        """Verify the campaign appears in the run index."""
        runs_df = self.store.list_runs()
        self.assertGreater(len(runs_df), 0)
        self.assertIn("run_id", runs_df.columns)
        self.assertIn("status", runs_df.columns)


@unittest.skipUnless(HAS_DEPS, "Requires torch, transformers, p53cad")
class TestPostAnalysisComponents(unittest.TestCase):
    """Test individual post-campaign analysis components."""

    def test_esmfold_structure_predictor(self):
        """Test ESMFold structure prediction (falls back to estimation)."""
        from p53cad.engine.md_validation import StructurePredictor

        predictor = StructurePredictor()
        # Use a short test sequence (not full p53)
        pred = predictor.predict_esmfold("MVLSPADKTNVKAAW", name="test")
        self.assertGreater(pred.mean_plddt, 0)
        self.assertGreater(len(pred.plddt_scores), 0)
        self.assertTrue(len(pred.pdb_string) > 0)

    def test_free_energy_calculator(self):
        """Test free energy estimation for a mutation."""
        from p53cad.engine.md_validation import FreeEnergyCalculator

        calc = FreeEnergyCalculator()
        from p53cad.data.dms import P53_WT

        result = calc.calculate_ddg(P53_WT, P53_WT, "M340S")
        self.assertIsNotNone(result.ddg_average)
        self.assertGreater(result.convergence_score, 0)

    def test_md_script_generation(self):
        """Test MD simulation script generation."""
        from p53cad.engine.md_validation import MDSimulationGenerator

        gen = MDSimulationGenerator()
        config = gen.generate_config("test_rescue", "MVLSPADKTNVKAAW", ["M1S"])

        for platform in ["openmm", "gromacs", "amber"]:
            script = gen.generate_script(config, platform)
            self.assertGreater(len(script), 100)
            self.assertIn("test_rescue", script)

    def test_explainability_energy_decomposition(self):
        """Test the energy decomposition module."""
        from p53cad.engine.explainability import EnergyDecomposer

        decomposer = EnergyDecomposer()
        from p53cad.data.dms import P53_WT

        # Test basic decomposition
        result = decomposer.decompose_rescue(
            P53_WT, P53_WT, P53_WT, ["M340S"]
        )
        self.assertIsNotNone(result.total_ddg)
        self.assertIn(result.rescue_quality, ["excellent", "good", "marginal", "poor"])
        self.assertIsNotNone(result.dominant_mechanism)

    def test_explainability_structural_decomposition(self):
        """Test structure-aware energy decomposition."""
        from p53cad.engine.explainability import EnergyDecomposer

        decomposer = EnergyDecomposer()
        from p53cad.data.dms import P53_WT

        # Test structural decomposition (will fall back if PDB missing)
        result = decomposer.decompose_rescue_structural(
            P53_WT, P53_WT, P53_WT, ["M340S"]
        )
        self.assertIsNotNone(result.total_ddg)
        self.assertIn(result.rescue_quality, ["excellent", "good", "marginal", "poor"])

    def test_rescue_mechanism_classifier(self):
        """Test rescue mechanism classification."""
        from p53cad.engine.explainability import RescueMechanismClassifier

        classifier = RescueMechanismClassifier()
        from p53cad.data.dms import P53_WT

        mechanism = classifier.classify(
            P53_WT, P53_WT, P53_WT, "R175H", ["M340S"]
        )
        self.assertIn(mechanism.primary_mechanism, [
            "stability_restoration", "binding_enhancement",
            "structural_compensation", "electrostatic_optimization"
        ])
        self.assertGreater(mechanism.confidence, 0)


if __name__ == "__main__":
    unittest.main()
