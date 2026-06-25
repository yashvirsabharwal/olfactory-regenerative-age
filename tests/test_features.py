from pathlib import Path
import sys
import unittest

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.features import (
    build_ora_feature_matrix,
    feature_kind_counts,
    merge_donor_feature_matrices,
    summarize_feature_family_stability,
)


class FeatureMatrixTests(unittest.TestCase):
    def test_build_feature_matrix_keeps_composition_only_by_default(self):
        features = pd.DataFrame(
            {
                "donor_id": ["d1"],
                "total_cells": [100],
                "prop__hbc": [0.2],
                "clr__mosn": [1.0],
                "ratio__lineage_fraction": [0.4],
                "technical": [5],
            }
        )

        matrix = build_ora_feature_matrix(features)

        self.assertEqual(matrix.columns.tolist(), ["donor_id", "prop__hbc", "clr__mosn", "ratio__lineage_fraction"])

    def test_build_feature_matrix_merges_module_scores(self):
        features = pd.DataFrame({"donor_id": ["d1", "d2"], "prop__hbc": [0.2, 0.4]})
        modules = pd.DataFrame(
            {
                "donor_id": ["d1", "d2"],
                "module_score__hbc_identity": [1.2, 1.4],
                "not_a_module": [9, 8],
            }
        )

        matrix = build_ora_feature_matrix(features, modules)
        counts = feature_kind_counts(matrix)

        self.assertEqual(matrix.columns.tolist(), ["donor_id", "prop__hbc", "module_score__hbc_identity"])
        self.assertEqual(
            counts,
            {"composition": 1, "module": 1, "scvi_global": 0, "scvi_cell_state": 0},
        )

    def test_merge_donor_feature_matrices_rejects_duplicate_features(self):
        left = pd.DataFrame({"donor_id": ["d1"], "prop__hbc": [0.2]})
        right = pd.DataFrame({"donor_id": ["d1"], "scvi_global_mean__dim01": [1.0]})
        merged = merge_donor_feature_matrices([left, right])
        self.assertEqual(merged.columns.tolist(), ["donor_id", "prop__hbc", "scvi_global_mean__dim01"])

        with self.assertRaises(ValueError):
            merge_donor_feature_matrices([left, left])

    def test_summarize_feature_family_stability_groups_hybrid_features(self):
        stability = pd.DataFrame(
            {
                "model": ["xgboost", "xgboost", "xgboost"],
                "feature": ["prop__hbc", "module_score__stress", "scvi_state_mean__hbc__dim01"],
                "abs_mean_importance": [0.2, 0.1, 0.4],
                "selection_fraction": [1.0, 0.5, 1.0],
            }
        )

        summary = summarize_feature_family_stability(stability)

        self.assertEqual(summary.iloc[0]["feature_family"], "scvi_cell_state")
        self.assertEqual(set(summary["feature_family"]), {"composition", "module", "scvi_cell_state"})


if __name__ == "__main__":
    unittest.main()
