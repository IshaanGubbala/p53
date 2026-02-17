import unittest

import numpy as np
import pandas as pd

from p53cad.results.visualization import (
    build_candidate_position_heatmap,
    build_loss_mesh,
    build_multi_metric_frame,
)


class ResultsVisualizationTests(unittest.TestCase):
    def test_candidate_heatmap_builder(self):
        matrix_df, freq_df = build_candidate_position_heatmap(
            [
                {"candidate_id": 1, "profile": "A", "score": 1.0, "mut_positions": [100, 150]},
                {"candidate_id": 2, "profile": "B", "score": 0.9, "mut_positions": [150, 200]},
            ],
            max_positions=10,
        )
        self.assertFalse(matrix_df.empty)
        self.assertFalse(freq_df.empty)
        self.assertIn("candidate_label", matrix_df.columns)

    def test_multi_metric_normalization(self):
        df = pd.DataFrame(
            {
                "step": [1, 2, 3],
                "score": [-0.2, 0.0, 0.4],
                "stability": [-0.2, -0.1, -0.05],
                "binding": [5.0, 6.5, 8.0],
                "identity": [90.0, 94.0, 96.0],
            }
        )
        out = build_multi_metric_frame(
            df,
            normalize=True,
            min_function=-0.2,
            min_stability=-0.2,
            min_binding=5.0,
            min_identity=90.0,
        )
        self.assertTrue(((out[["score", "stability", "binding", "identity"]] >= 0.0) & (out[["score", "stability", "binding", "identity"]] <= 1.0)).all().all())

    def test_loss_mesh_builder(self):
        rng = np.random.default_rng(0)
        df = pd.DataFrame(
            {
                "lx": rng.normal(size=40),
                "ly": rng.normal(size=40),
                "loss_total": rng.normal(loc=1.0, scale=0.2, size=40),
            }
        )
        mesh = build_loss_mesh(df, grid_size=20)
        self.assertIsNotNone(mesh)
        gx, gy, gz = mesh
        self.assertEqual(gx.shape, gy.shape)
        self.assertEqual(gx.shape, gz.shape)


if __name__ == "__main__":
    unittest.main()
