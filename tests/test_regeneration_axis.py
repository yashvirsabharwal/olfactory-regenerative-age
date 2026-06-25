import unittest

import numpy as np
import pandas as pd

from ora.regeneration_axis import (
    build_regeneration_axis_feature_map,
    build_regeneration_axis_theme_summary,
    build_regeneration_feature_resource_map,
    classify_regeneration_feature,
)


class RegenerationAxisTests(unittest.TestCase):
    def test_classify_regeneration_feature_maps_core_themes(self):
        self.assertEqual(
            classify_regeneration_feature("clr__quiescent_hbc").primary_theme,
            "basal_quiescence",
        )
        self.assertEqual(
            classify_regeneration_feature("module_score__immature_neuron").primary_theme,
            "immature_osn",
        )
        self.assertEqual(
            classify_regeneration_feature("prop__goblet").primary_theme,
            "respiratory_metaplasia_ciliated_goblet",
        )
        self.assertEqual(
            classify_regeneration_feature("clr__macrophage").primary_theme,
            "immune_inflammatory",
        )

    def test_resource_map_covers_all_numeric_features(self):
        matrix = pd.DataFrame(
            {
                "donor_id": ["d1", "d2"],
                "clr__quiescent_hbc": [0.1, 0.2],
                "prop__goblet": [0.3, 0.4],
                "note": ["a", "b"],
            }
        )

        resource = build_regeneration_feature_resource_map(matrix)

        self.assertEqual(set(resource["feature"]), {"clr__quiescent_hbc", "prop__goblet"})
        self.assertTrue(resource["primary_theme"].notna().all())
        self.assertTrue(resource["evidence_citations"].str.len().gt(0).all())

    def test_feature_map_joins_age_importance_and_cross_tissue_evidence(self):
        donors = [f"d{i}" for i in range(10)]
        matrix = pd.DataFrame(
            {
                "donor_id": donors,
                "clr__quiescent_hbc": np.linspace(0.0, 1.0, 10),
                "module_score__immature_neuron": np.linspace(1.0, 0.0, 10),
                "prop__goblet": np.linspace(0.2, 0.8, 10),
            }
        )
        manifest = pd.DataFrame(
            {
                "donor_id": donors,
                "age": np.linspace(25, 85, 10),
                "usable_for_ora_training": [True] * 10,
            }
        )
        stability = pd.DataFrame(
            {
                "model": ["xgboost", "catboost"],
                "feature": ["clr__quiescent_hbc", "prop__goblet"],
                "abs_mean_importance": [2.0, 1.0],
                "selection_fraction": [1.0, 0.5],
            }
        )
        cross_tissue = pd.DataFrame(
            {
                "feature": ["prop__goblet"],
                "specificity_class": ["airway_nasal_shared"],
                "classification_confidence": ["high"],
                "external_age_effect_status": ["pending"],
            }
        )

        feature_map = build_regeneration_axis_feature_map(
            feature_matrix=matrix,
            manifest=manifest,
            feature_stability=stability,
            cross_tissue_classification=cross_tissue,
        )
        summary = build_regeneration_axis_theme_summary(feature_map)

        hbc = feature_map[feature_map["feature"].eq("clr__quiescent_hbc")].iloc[0]
        goblet = feature_map[feature_map["feature"].eq("prop__goblet")].iloc[0]
        self.assertEqual(hbc["observed_vs_expected"], "aligned")
        self.assertEqual(float(hbc["max_abs_importance"]), 2.0)
        self.assertEqual(goblet["specificity_class"], "airway_nasal_shared")
        self.assertEqual(int(summary["n_features"].sum()), 3)


if __name__ == "__main__":
    unittest.main()
