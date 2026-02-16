"""
Tests for p53cad.engine.physics_validation

Unit tests run unconditionally; integration tests requiring openmm/mdtraj/transformers
are guarded by @unittest.skipUnless.
"""

import json
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

from p53cad.engine.physics_validation import (
    DDGResult,
    DNABindingInterfaceResult,
    EnergyMinimizationResult,
    ESMFoldResult,
    MDStabilityResult,
    PhysicsValidationPipeline,
    PhysicsValidationReport,
    PhysicsValidationResult,
    compute_physics_score,
    _seq_hash,
)


class TestSeqHash(unittest.TestCase):
    def test_deterministic(self):
        h1 = _seq_hash("MEEPQ")
        h2 = _seq_hash("MEEPQ")
        self.assertEqual(h1, h2)

    def test_different_seqs(self):
        self.assertNotEqual(_seq_hash("MEEPQ"), _seq_hash("MEEPX"))

    def test_length(self):
        self.assertEqual(len(_seq_hash("MEEPQ")), 16)


class TestDataclassSerialization(unittest.TestCase):
    """Dataclass to_dict round-trip — pdb_string must be excluded."""

    def test_esmfold_excludes_pdb(self):
        r = ESMFoldResult(
            pdb_string="HEADER ...\nEND",
            plddt_scores=[80.0, 85.0],
            mean_plddt=82.5,
            dbd_plddt=84.0,
        )
        d = r.to_dict()
        self.assertNotIn("pdb_string", d)
        self.assertAlmostEqual(d["mean_plddt"], 82.5)
        self.assertAlmostEqual(d["dbd_plddt"], 84.0)

    def test_energy_round_trip(self):
        r = EnergyMinimizationResult(potential_energy_kcal=-5000.0, n_atoms=6000)
        d = r.to_dict()
        self.assertEqual(d["potential_energy_kcal"], -5000.0)
        self.assertEqual(d["n_atoms"], 6000)
        # Should be JSON-serializable
        json.dumps(d)

    def test_ddg_round_trip(self):
        r = DDGResult(ddg_vs_wt_kcal=-8.5, ddg_vs_cancer_kcal=-15.0, interpretation="stabilizing")
        d = r.to_dict()
        self.assertEqual(d["interpretation"], "stabilizing")
        json.dumps(d)

    def test_md_round_trip(self):
        r = MDStabilityResult(
            rmsd_mean=1.5, rmsd_std=0.3, rmsd_final=1.8,
            rmsf_per_residue=[0.5] * 393,
            rmsf_dbd_mean=0.6, ss_content={"helix": 0.3, "sheet": 0.2, "coil": 0.5},
            rg_mean=1.8, stability_score=80.0, stability_verdict="stable",
        )
        d = r.to_dict()
        self.assertEqual(d["stability_verdict"], "stable")
        self.assertEqual(len(d["rmsf_per_residue"]), 393)
        json.dumps(d)

    def test_dna_binding_round_trip(self):
        r = DNABindingInterfaceResult(
            dbd_rmsd=1.2,
            per_contact_rmsd={"K120": 0.5, "R248": 0.8},
            interface_preservation_score=0.85,
        )
        d = r.to_dict()
        self.assertEqual(d["per_contact_rmsd"]["R248"], 0.8)
        json.dumps(d)

    def test_physics_result_composite(self):
        r = PhysicsValidationResult(
            candidate_rank=1,
            target_label="R175H+gene_therapy",
            sequence_hash="abc123",
            esmfold=ESMFoldResult("pdb", [80], 80, 82),
            ddg=DDGResult(-5, -10, "stabilizing"),
            overall_physics_score=72.5,
            verdict="PROMISING",
        )
        d = r.to_dict()
        self.assertIn("esmfold", d)
        self.assertNotIn("pdb_string", d.get("esmfold", {}))
        self.assertEqual(d["verdict"], "PROMISING")
        json.dumps(d)

    def test_report_to_dict(self):
        report = PhysicsValidationReport(
            run_id="test_run",
            n_candidates=2,
            n_tier2=0,
            wt_energy_kcal=-5000.0,
            candidates=[{"rank": 1}, {"rank": 2}],
        )
        d = report.to_dict()
        self.assertEqual(d["run_id"], "test_run")
        self.assertEqual(len(d["candidates"]), 2)
        json.dumps(d)


class TestCompositeScoring(unittest.TestCase):
    """Test the 0-100 composite physics score logic."""

    def test_all_none_gives_neutral(self):
        score, verdict = compute_physics_score()
        # All neutral fallbacks: 15+12.5+12.5+10 = 50
        self.assertAlmostEqual(score, 50.0, places=0)
        self.assertIn(verdict, ("PROMISING", "UNCERTAIN"))

    def test_perfect_scores(self):
        ef = ESMFoldResult("", [], 90, 90)
        ddg = DDGResult(-15, -20, "stabilizing")
        md = MDStabilityResult(0.5, 0.1, 0.5, [], 0.3, {}, 1.5, 100, "stable")
        dna = DNABindingInterfaceResult(0.5, {}, 1.0)
        score, verdict = compute_physics_score(ef, ddg, md, dna)
        self.assertGreater(score, 90)
        self.assertEqual(verdict, "STRONG")

    def test_terrible_scores(self):
        ef = ESMFoldResult("", [], 25, 25)
        ddg = DDGResult(20, 15, "destabilizing")
        md = MDStabilityResult(5.0, 1.0, 6.0, [], 3.0, {}, 2.5, 0, "unstable")
        dna = DNABindingInterfaceResult(8.0, {}, 0.0)
        score, verdict = compute_physics_score(ef, ddg, md, dna)
        self.assertLess(score, 20)
        self.assertEqual(verdict, "CONCERNING")

    def test_partial_data(self):
        """Only ESMFold available, rest missing."""
        ef = ESMFoldResult("", [], 80, 80)
        score, verdict = compute_physics_score(esmfold=ef)
        # ESMFold: (80-30)/60 * 30 = 25, rest neutral: 12.5+12.5+10 = 35
        # total = 60/100 * 100 = 60
        self.assertGreater(score, 50)
        self.assertIn(verdict, ("PROMISING", "STRONG"))

    def test_verdict_boundaries(self):
        """Test verdict string thresholds."""
        # Create results that target specific score ranges
        # STRONG >= 75
        ef_strong = ESMFoldResult("", [], 90, 90)
        ddg_strong = DDGResult(-10, -15, "stabilizing")
        dna_strong = DNABindingInterfaceResult(0.5, {}, 0.95)
        s, v = compute_physics_score(esmfold=ef_strong, ddg=ddg_strong, dna_binding=dna_strong)
        self.assertEqual(v, "STRONG")

        # CONCERNING < 35
        ef_bad = ESMFoldResult("", [], 30, 30)
        ddg_bad = DDGResult(15, 10, "destabilizing")
        dna_bad = DNABindingInterfaceResult(6.0, {}, 0.0)
        s, v = compute_physics_score(esmfold=ef_bad, ddg=ddg_bad, dna_binding=dna_bad)
        self.assertEqual(v, "CONCERNING")


class TestPipelineWiring(unittest.TestCase):
    """Test pipeline orchestration with all steps skipped (verifies wiring only)."""

    def test_all_skipped(self):
        import pandas as pd

        pipeline = PhysicsValidationPipeline(device="cpu")
        top_df = pd.DataFrame({
            "sequence": ["MEEPQ" * 78 + "MSD"],  # 393 residues
            "target_label": ["R175H+gene_therapy"],
        })

        report = pipeline.validate_campaign(
            run_id="test_skip",
            top_df=top_df,
            output_dir=Path("/tmp/test_physics_skip"),
            wt_sequence="MEEPQ" * 78 + "MSD",
            skip_esmfold=True,
            skip_energy=True,
            skip_md=True,
            skip_dna=True,
        )

        self.assertEqual(report.n_candidates, 1)
        self.assertEqual(report.n_tier2, 0)
        self.assertIsNone(report.wt_energy_kcal)
        self.assertEqual(len(report.candidates), 1)
        c = report.candidates[0]
        self.assertEqual(c["target_label"], "R175H+gene_therapy")
        # All steps skipped -> neutral composite score
        self.assertAlmostEqual(c["overall_physics_score"], 50.0, places=0)


# ---------------------------------------------------------------------------
# Integration tests (require real dependencies)
# ---------------------------------------------------------------------------

def _has_openmm():
    try:
        import openmm
        return True
    except ImportError:
        return False


def _has_mdtraj():
    try:
        import mdtraj
        return True
    except ImportError:
        return False


def _has_esmfold():
    try:
        from transformers import EsmForProteinFolding
        return True
    except ImportError:
        return False


@unittest.skipUnless(_has_openmm(), "OpenMM not installed")
class TestOpenMMIntegration(unittest.TestCase):
    """Real OpenMM energy minimization on a tiny PDB snippet."""

    def test_minimize_small_pdb(self):
        from p53cad.engine.physics_validation import OpenMMEnergyCalculator

        # Minimal alanine dipeptide-style PDB
        pdb = """HEADER    TEST
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1       1.246   2.390   0.000  1.00  0.00           O
ATOM      5  CB  ALA A   1       1.986  -0.764   1.193  1.00  0.00           C
ATOM      6  N   ALA A   2       3.323   1.510   0.000  1.00  0.00           N
ATOM      7  CA  ALA A   2       3.979   2.812   0.000  1.00  0.00           C
ATOM      8  C   ALA A   2       5.487   2.686   0.000  1.00  0.00           C
ATOM      9  O   ALA A   2       6.104   1.610   0.000  1.00  0.00           O
ATOM     10  CB  ALA A   2       3.543   3.635  -1.213  1.00  0.00           C
END
"""
        calc = OpenMMEnergyCalculator()
        result = calc.minimize_and_get_energy(pdb, max_iterations=50)
        self.assertIsInstance(result.potential_energy_kcal, float)
        self.assertGreater(result.n_atoms, 0)

    def test_ddg_computation(self):
        from p53cad.engine.physics_validation import OpenMMEnergyCalculator

        calc = OpenMMEnergyCalculator()
        ddg = calc.compute_ddg(e_wt=-5000, e_cancer=-4900, e_rescue=-4950)
        self.assertEqual(ddg.ddg_vs_wt_kcal, 50.0)  # rescue - wt
        self.assertEqual(ddg.ddg_vs_cancer_kcal, -50.0)  # rescue - cancer
        self.assertEqual(ddg.interpretation, "destabilizing")  # >+5


@unittest.skipUnless(_has_openmm() and _has_mdtraj(), "OpenMM + mdtraj not installed")
class TestMDIntegration(unittest.TestCase):
    """Real short MD simulation on a tiny system."""

    def test_very_short_md(self):
        from p53cad.engine.physics_validation import MDStabilityChecker

        # Use a small alanine chain for speed
        pdb = "HEADER    TEST\n"
        for i in range(5):
            x = i * 3.8
            pdb += (
                f"ATOM  {i+1:5d}  CA  ALA A{i+1:4d}    "
                f"{x:8.3f}   0.000   0.000  1.00 50.00           C\n"
            )
        pdb += "END\n"

        checker = MDStabilityChecker()
        # 0.001 ns = 1 ps, extremely fast
        result = checker.run_stability_check(pdb, simulation_ns=0.001, name="tiny_test")
        self.assertIsInstance(result.rmsd_mean, float)
        self.assertIn(result.stability_verdict, ("stable", "metastable", "unstable"))


class TestRuntimeProbes(unittest.TestCase):
    """Verify physics validation probes exist in runtime capabilities."""

    def test_probes_present(self):
        from p53cad.core.runtime import get_runtime_capabilities

        caps = get_runtime_capabilities()
        self.assertIn("pdbfixer_installed", caps)
        self.assertIn("mdtraj_installed", caps)
        self.assertIn("esmfold_local_available", caps)
        self.assertIn("physics_validation_ready", caps)
        # All should be boolean
        self.assertIsInstance(caps["physics_validation_ready"], bool)


if __name__ == "__main__":
    unittest.main()
