from pathlib import Path
import sys
import unittest

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.model_compare import compare_feature_set_deltas, rank_feature_set_summaries


class ModelCompareTests(unittest.TestCase):
    def test_rank_feature_set_summaries_marks_best(self):
        base = _summary({"ridge": 15.0, "random_forest": 14.2})
        augmented = _summary({"ridge": 14.8, "random_forest": 14.0})

        ranked = rank_feature_set_summaries({"composition": base, "composition_plus_modules": augmented})

        self.assertEqual(ranked.iloc[0]["feature_set"], "composition_plus_modules")
        self.assertEqual(ranked.iloc[0]["model"], "random_forest")
        self.assertTrue(bool(ranked.iloc[0]["is_best_overall"]))
        self.assertEqual(set(ranked["mae_rank_within_feature_set"]), {1, 2})
        self.assertIn("backend_package", ranked.columns)

    def test_compare_feature_set_deltas_uses_augmented_minus_base(self):
        base = _summary({"ridge": 15.0, "random_forest": 14.2})
        augmented = _summary({"ridge": 14.8, "random_forest": 14.4})

        deltas = compare_feature_set_deltas(base, augmented)

        ridge = deltas[deltas["model"].eq("ridge")].iloc[0]
        forest = deltas[deltas["model"].eq("random_forest")].iloc[0]
        self.assertAlmostEqual(ridge["delta_mae_mean"], -0.2)
        self.assertTrue(bool(ridge["mae_improved"]))
        self.assertAlmostEqual(forest["delta_mae_mean"], 0.2)
        self.assertFalse(bool(forest["mae_improved"]))


def _summary(mae_by_model: dict[str, float]) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "model": model,
                "mae_mean": mae,
                "rmse_mean": mae + 3.0,
                "r2_mean": 0.1,
                "spearman_r_mean": 0.3,
                "backend_package": "scikit-learn",
                "backend_version": "1.5.0",
                "fallback_used": False,
            }
            for model, mae in mae_by_model.items()
        ]
    )


if __name__ == "__main__":
    unittest.main()
