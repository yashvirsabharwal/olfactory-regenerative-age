import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import anndata as ad
import numpy as np
import pandas as pd

from ora.scvi_donor import build_scvi_donor_embedding_features, summarize_scvi_state_importance


class ScviDonorTests(unittest.TestCase):
    def test_build_scvi_donor_embedding_features_aggregates_global_and_state_means(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "toy.h5ad"
            adata = ad.AnnData(X=np.ones((6, 3)))
            adata.obs["donor_id"] = ["d1", "d1", "d1", "d2", "d2", "d2"]
            adata.obs["fine_celltype"] = ["state A", "state A", "state B", "state A", "state B", "state B"]
            adata.obsm["X_scvi"] = np.array(
                [
                    [1.0, 2.0],
                    [3.0, 4.0],
                    [5.0, 6.0],
                    [2.0, 4.0],
                    [4.0, 8.0],
                    [6.0, 12.0],
                ]
            )
            adata.write_h5ad(path)

            features, qc = build_scvi_donor_embedding_features(
                path,
                top_cell_states=2,
                min_cells_per_state=2,
                chunk_size=2,
            )

        d1 = features.set_index("donor_id").loc["d1"]
        d2 = features.set_index("donor_id").loc["d2"]
        self.assertAlmostEqual(d1["scvi_global_mean__dim01"], 3.0)
        self.assertAlmostEqual(d2["scvi_global_mean__dim02"], 8.0)
        self.assertAlmostEqual(d1["scvi_state_mean__state_a__dim01"], 2.0)
        self.assertTrue(np.isnan(d2["scvi_state_mean__state_a__dim01"]))
        self.assertAlmostEqual(d2["scvi_state_mean__state_b__dim02"], 10.0)
        self.assertIn("missing_fraction_at_min", qc.columns)

    def test_summarize_scvi_state_importance_groups_global_and_cell_state_features(self):
        feature_stability = np.array(
            [
                ["catboost", "scvi_global_mean__dim01", 0.5, 1.0],
                ["catboost", "scvi_state_mean__fully_mature_mosn__dim02", 0.3, 0.8],
                ["catboost", "unrelated_feature", 1.0, 1.0],
            ],
            dtype=object,
        )
        summary = summarize_scvi_state_importance(
            pd.DataFrame(
                feature_stability,
                columns=["model", "feature", "abs_mean_importance", "selection_fraction"],
            )
        )

        self.assertEqual(summary.shape[0], 2)
        self.assertIn("all_cells", summary["cell_state"].tolist())
        self.assertIn("fully_mature_mosn", summary["cell_state"].tolist())
        global_row = summary[summary["cell_state"].eq("all_cells")].iloc[0]
        self.assertEqual(global_row["state_rank_within_model"], 1)


if __name__ == "__main__":
    unittest.main()
