import unittest

import numpy as np
import pandas as pd

from p53cad.results.schema import select_presentation_shortlist


class ShortlistIntegrityTests(unittest.TestCase):
    def _base_rows(self):
        return [
            {
                "candidate_uid": "valid_1",
                "scenario_id": "R175H__gene_therapy",
                "target_label": "R175H",
                "targets_json": "[\"R175H\"]",
                "delivery_method": "gene_therapy",
                "profile": "Balanced",
                "score": 1.2,
                "clinical_score": 80,
                "uncertainty": 0.1,
                "ood_distance": 0.8,
                "n_mutations": 2,
                "mutations_json": "[\"R175H\",\"N239Y\"]",
                "sequence": "SEQ_VALID_1",
                "pareto_rank": 1,
            },
            {
                "candidate_uid": "invalid_missing_target",
                "scenario_id": "R175H__mrna_therapy",
                "target_label": "R175H",
                "targets_json": "[\"R175H\"]",
                "delivery_method": "mrna_therapy",
                "profile": "Stability-First",
                "score": 1.3,
                "clinical_score": 82,
                "uncertainty": 0.11,
                "ood_distance": 0.7,
                "n_mutations": 1,
                "mutations_json": "[\"N239Y\"]",
                "sequence": "SEQ_INVALID_NO_TARGET",
                "pareto_rank": 1,
            },
            {
                "candidate_uid": "invalid_zero_mut",
                "scenario_id": "R248Q__protein_therapy",
                "target_label": "R248Q",
                "targets_json": "[\"R248Q\"]",
                "delivery_method": "protein_therapy",
                "profile": "Conservative",
                "score": 1.4,
                "clinical_score": 84,
                "uncertainty": 0.08,
                "ood_distance": 0.6,
                "n_mutations": 0,
                "mutations_json": "[]",
                "sequence": "SEQ_INVALID_ZERO",
                "pareto_rank": 2,
            },
            {
                "candidate_uid": "valid_2",
                "scenario_id": "R248Q__protein_therapy",
                "target_label": "R248Q",
                "targets_json": "[\"R248Q\"]",
                "delivery_method": "protein_therapy",
                "profile": "Experimental",
                "score": 1.1,
                "clinical_score": 79,
                "uncertainty": 0.2,
                "ood_distance": 1.0,
                "n_mutations": 2,
                "mutations_json": "[\"R248Q\",\"S241H\"]",
                "sequence": "SEQ_VALID_2",
                "pareto_rank": 2,
            },
        ]

    def test_strict_retain_add_excludes_invalid_candidates(self):
        df = pd.DataFrame(self._base_rows())
        out = select_presentation_shortlist(df, top_n=10)
        self.assertFalse(out.empty)
        self.assertTrue((out["n_mutations"] > 0).all())
        self.assertNotIn("invalid_missing_target", set(out["candidate_uid"]))
        self.assertNotIn("invalid_zero_mut", set(out["candidate_uid"]))

    def test_dedup_keeps_best_representative(self):
        rows = self._base_rows()
        rows.append(
            {
                "candidate_uid": "dup_lower",
                "scenario_id": "R175H__gene_therapy",
                "target_label": "R175H",
                "targets_json": "[\"R175H\"]",
                "delivery_method": "gene_therapy",
                "profile": "Balanced",
                "score": 1.0,
                "clinical_score": 70,
                "uncertainty": 0.2,
                "ood_distance": 1.1,
                "n_mutations": 2,
                "mutations_json": "[\"R175H\",\"N239Y\"]",
                "sequence": "SEQ_VALID_1",
                "pareto_rank": 2,
            }
        )
        out = select_presentation_shortlist(pd.DataFrame(rows), top_n=10)
        uid_set = set(out["candidate_uid"])
        self.assertIn("valid_1", uid_set)
        self.assertNotIn("dup_lower", uid_set)

    def test_selection_metadata_is_finite(self):
        df = pd.DataFrame(self._base_rows())
        out = select_presentation_shortlist(df, top_n=10)
        for col in ["selection_score", "diversity_novelty", "max_mutation_overlap"]:
            self.assertIn(col, out.columns)
            vals = pd.to_numeric(out[col], errors="coerce")
            self.assertTrue(np.isfinite(vals).all())
        self.assertTrue(out["selection_reason"].astype(str).str.len().gt(0).all())

    def test_target_and_delivery_diversity(self):
        """Verify shortlist spreads across targets and delivery methods."""
        hotspots = ["R175H", "R248Q", "R248W", "R273H", "G245S", "R249S", "R282W"]
        deliveries = ["gene_therapy", "mrna_therapy", "protein_therapy"]
        rows = []
        uid = 0
        for h in hotspots:
            for d in deliveries:
                # 3 candidates per hotspot-delivery combo, decreasing score
                for k in range(3):
                    uid += 1
                    rows.append({
                        "candidate_uid": f"cand_{uid}",
                        "scenario_id": f"{h}__{d}",
                        "target_label": h,
                        "targets_json": f'["{h}"]',
                        "delivery_method": d,
                        "profile": ["Balanced", "Stability-First", "Conservative"][k],
                        "score": 1.5 - uid * 0.005,
                        "clinical_score": 80 - uid * 0.5,
                        "uncertainty": 0.05 + uid * 0.002,
                        "ood_distance": 0.3 + uid * 0.01,
                        "n_mutations": 2,
                        "mutations_json": f'["{h}","A{100+uid}G"]',
                        "sequence": f"SEQ_{uid}",
                        "pareto_rank": 1 + k,
                    })
        df = pd.DataFrame(rows)
        out = select_presentation_shortlist(df, top_n=30, max_per_target=4)
        self.assertGreater(len(out), 0)

        # At least 5 distinct targets in top-15
        top15 = out.head(15)
        self.assertGreaterEqual(top15["target_label"].nunique(), 5,
            f"Expected >=5 targets in top-15, got {top15['target_label'].nunique()}")

        # All 3 delivery methods in top-30
        delivery_set = set(out["delivery_method"].str.lower())
        for d in deliveries:
            self.assertIn(d, delivery_set,
                f"Delivery method '{d}' missing from top-30")

    def test_base_mutation_cap(self):
        """Verify no single hotspot dominates via combo targets."""
        hotspots = ["R175H", "R248Q", "R248W", "R273H", "G245S", "R249S"]
        rows = []
        uid = 0
        # Create many combos containing R249S
        for h in hotspots:
            for other in hotspots:
                if h >= other:
                    continue
                uid += 1
                combo = "+".join(sorted([h, other]))
                muts = sorted([h, other, f"A{200+uid}G"])
                rows.append({
                    "candidate_uid": f"combo_{uid}",
                    "scenario_id": f"{combo}__gene_therapy",
                    "target_label": combo,
                    "targets_json": f'["{h}","{other}"]',
                    "delivery_method": "gene_therapy",
                    "profile": "Balanced",
                    "score": 1.5 - uid * 0.01,
                    "clinical_score": 80,
                    "uncertainty": 0.1,
                    "ood_distance": 0.5,
                    "n_mutations": 3,
                    "mutations_json": str(muts).replace("'", '"'),
                    "sequence": f"SEQ_COMBO_{uid}",
                    "pareto_rank": 1,
                })
        df = pd.DataFrame(rows)
        out = select_presentation_shortlist(df, top_n=15, max_per_base_mutation=4)

        # Count base mutation appearances
        from collections import Counter
        base_counts = Counter()
        for _, row in out.iterrows():
            label = str(row["target_label"])
            for m in label.split("+"):
                base_counts[m.strip()] += 1

        for mut, cnt in base_counts.items():
            self.assertLessEqual(cnt, 4,
                f"Base mutation {mut} appeared {cnt} times, exceeds cap of 4")


if __name__ == "__main__":
    unittest.main()
