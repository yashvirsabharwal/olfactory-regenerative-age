from pathlib import Path
import sys
import unittest

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.manuscript import build_model_card


class ManuscriptTests(unittest.TestCase):
    def test_build_model_card_combines_benchmark_context(self):
        comparison = pd.DataFrame(
            {
                "feature_set": ["composition_plus_modules", "composition"],
                "model": ["catboost", "xgboost"],
                "n": [187, 187],
                "repeats": [20, 20],
                "mae_mean": [14.08, 14.15],
                "mae_ci_low": [13.6, 13.7],
                "mae_ci_high": [14.5, 14.7],
                "spearman_r_mean": [0.34, 0.35],
                "is_best_overall": [True, False],
            }
        )
        calibration = pd.DataFrame(
            {"model": ["catboost"], "calibration_slope_ora_on_age": [0.15]}
        )
        permutation = pd.DataFrame({"model": ["catboost"], "empirical_p_mae": [0.02]})

        card = build_model_card(
            feature_set_comparison=comparison,
            calibration=calibration,
            permutation=permutation,
        )

        best = card[card["model"].eq("catboost")].iloc[0]
        self.assertEqual(best["role"], "preferred_benchmark")
        self.assertEqual(best["permutation_p_mae"], 0.02)
        self.assertIn("module gain is modest", best["limitations"])


if __name__ == "__main__":
    unittest.main()
