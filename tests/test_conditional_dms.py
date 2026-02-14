"""Tests for conditional DMS (rescue) scoring.

Validates that compute_conditional_rescue_scores returns finite log-probs
with correct shapes for masked positions in cancer-context sequences.

NOTE: Requires PyTorch and ESM-2 model.
"""

import os
import unittest

os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

try:
    import torch
    from p53cad.core.runtime import bootstrap_runtime

    bootstrap_runtime()
    from p53cad.data.dms import P53_WT, apply_mutation
    from p53cad.engine.latent import ManifoldEmbedder
    from p53cad.engine.oracle import compute_conditional_rescue_scores

    HAS_DEPS = True
except ImportError:
    HAS_DEPS = False


STANDARD_AAS = "ACDEFGHIKLMNPQRSTVWY"


@unittest.skipUnless(HAS_DEPS, "Requires torch, transformers, p53cad")
class TestConditionalRescueScores(unittest.TestCase):
    """Test compute_conditional_rescue_scores function."""

    @classmethod
    def setUpClass(cls):
        cls.embedder = ManifoldEmbedder()

    def test_returns_dict_for_each_position(self):
        """Result has an entry for each requested position."""
        cancer_seq = apply_mutation(P53_WT, "R175H")
        positions = [200, 220, 240]
        result = compute_conditional_rescue_scores(
            self.embedder, cancer_seq, positions
        )
        self.assertIsInstance(result, dict)
        for pos in positions:
            self.assertIn(pos, result)

    def test_all_20_amino_acids_present(self):
        """Each position has log-prob entries for all 20 standard amino acids."""
        cancer_seq = apply_mutation(P53_WT, "R248Q")
        positions = [180]
        result = compute_conditional_rescue_scores(
            self.embedder, cancer_seq, positions
        )
        aa_dict = result[180]
        for aa in STANDARD_AAS:
            self.assertIn(aa, aa_dict, f"Missing amino acid {aa} at position 180")

    def test_values_are_finite(self):
        """All returned log-probs should be finite numbers."""
        import math

        cancer_seq = apply_mutation(P53_WT, "G245S")
        positions = [150, 250]
        result = compute_conditional_rescue_scores(
            self.embedder, cancer_seq, positions
        )
        for pos, aa_scores in result.items():
            for aa, logp in aa_scores.items():
                self.assertFalse(
                    math.isnan(logp), f"NaN log-prob at pos={pos}, aa={aa}"
                )
                self.assertFalse(
                    math.isinf(logp), f"Inf log-prob at pos={pos}, aa={aa}"
                )

    def test_log_probs_are_negative(self):
        """Log-probabilities should be non-positive (log of [0,1])."""
        cancer_seq = apply_mutation(P53_WT, "R175H")
        positions = [200]
        result = compute_conditional_rescue_scores(
            self.embedder, cancer_seq, positions
        )
        for aa, logp in result[200].items():
            self.assertLessEqual(logp, 0.0, f"Positive log-prob at aa={aa}: {logp}")

    def test_empty_positions_returns_empty(self):
        """Empty position list returns empty dict."""
        cancer_seq = apply_mutation(P53_WT, "R175H")
        result = compute_conditional_rescue_scores(
            self.embedder, cancer_seq, []
        )
        self.assertEqual(len(result), 0)


if __name__ == "__main__":
    unittest.main()
