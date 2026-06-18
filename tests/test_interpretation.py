from pathlib import Path
import sys
import unittest

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.interpretation import build_feature_interpretation, classify_feature_theme


class FeatureInterpretationTests(unittest.TestCase):
    def test_classify_feature_theme_maps_core_biology(self):
        self.assertEqual(classify_feature_theme("clr__quiescent_hbc"), "regenerative/progenitor epithelium")
        self.assertEqual(classify_feature_theme("module_score__immature_neuron"), "neuronal lineage maturation")
        self.assertEqual(classify_feature_theme("clr__macrophage"), "immune/inflammatory compartment")
        self.assertEqual(classify_feature_theme("prop__mucous_gland"), "supporting/secretory epithelium")

    def test_build_feature_interpretation_aggregates_model_support(self):
        stability = pd.DataFrame(
            {
                "model": ["catboost", "xgboost", "catboost", "null_model"],
                "feature": [
                    "clr__quiescent_hbc",
                    "clr__quiescent_hbc",
                    "module_score__immature_neuron",
                    "clr__macrophage",
                ],
                "selection_fraction": [1.0, 0.8, 1.0, 1.0],
                "abs_mean_importance": [2.0, 1.0, 3.0, 10.0],
            }
        )
        associations = pd.DataFrame(
            {
                "feature": ["clr__quiescent_hbc"],
                "direction": ["positive"],
                "beta_per_10_years": [0.12],
                "fdr": [0.01],
            }
        )

        result = build_feature_interpretation(stability, associations, top_per_model=3, top_n=5)

        hbc = result[result["feature"].eq("clr__quiescent_hbc")].iloc[0]
        self.assertEqual(hbc["n_supporting_models"], 2)
        self.assertEqual(hbc["age_direction"], "positive")
        self.assertIn("catboost", hbc["supporting_models"])
        self.assertNotIn("null_model", result["supporting_models"].str.cat(sep=","))


if __name__ == "__main__":
    unittest.main()
