import unittest

import numpy as np
import pandas as pd

from ora.cross_tissue import (
    build_cross_tissue_candidate_matrix,
    build_cross_tissue_specificity_summary,
    build_ora_cross_tissue_feature_classification,
    classify_feature,
)


class CrossTissueSpecificityTests(unittest.TestCase):
    def test_candidate_matrix_selects_airway_lung_specificity_resources_and_placeholders(self):
        config = {
            "public_data_exhaustion": {
                "candidates": [
                    {
                        "accession_or_dataset": "LungMAP",
                        "source_url": "https://example.org/lungmap",
                        "tissue": "Lung",
                        "assay": "snRNA-seq",
                        "species": "human",
                        "age_availability": "yes",
                        "inclusion_decision": "include for healthy aging airway/lung comparator",
                    },
                    {
                        "accession_or_dataset": "Mouse OE",
                        "tissue": "Main olfactory epithelium",
                        "assay": "scRNA-seq",
                        "species": "mouse",
                        "inclusion_decision": "exclude",
                    },
                ]
            }
        }

        matrix = build_cross_tissue_candidate_matrix(config)

        self.assertIn("lungmap", set(matrix["dataset_id"]))
        self.assertIn("skin_epithelium_query_required", set(matrix["dataset_id"]))
        self.assertNotIn("mouse_oe", set(matrix["dataset_id"]))

    def test_classify_feature_rules_cover_main_specificity_classes(self):
        self.assertEqual(classify_feature("clr__early_iosn").specificity_class, "olfactory_specific")
        self.assertEqual(classify_feature("prop__goblet").specificity_class, "airway_nasal_shared")
        self.assertEqual(
            classify_feature("module_score__hbc_activation_injury").specificity_class,
            "pan_epithelial_regenerative",
        )
        self.assertEqual(classify_feature("clr__macrophage").specificity_class, "immune_inflammatory_shared")

    def test_feature_classification_joins_age_and_importance_evidence(self):
        donors = [f"d{i}" for i in range(8)]
        features = pd.DataFrame(
            {
                "donor_id": donors,
                "clr__early_iosn": np.linspace(1, 0, 8),
                "prop__goblet": np.linspace(0, 1, 8),
                "module_score__hbc_activation_injury": np.linspace(0.2, 0.8, 8),
            }
        )
        manifest = pd.DataFrame(
            {
                "donor_id": donors,
                "age": np.linspace(30, 80, 8),
                "usable_for_ora_training": [True] * 8,
            }
        )
        stability = pd.DataFrame(
            {
                "model": ["ridge", "ridge"],
                "feature": ["clr__early_iosn", "prop__goblet"],
                "abs_mean_importance": [2.0, 1.0],
                "selection_fraction": [1.0, 0.5],
            }
        )
        comparators = pd.DataFrame(
            {
                "dataset_id": ["lungmap"],
                "tissue_class": ["lung"],
                "local_status": ["candidate_selected_external_effects_pending"],
            }
        )

        classification = build_ora_cross_tissue_feature_classification(
            feature_matrix=features,
            manifest=manifest,
            feature_stability=stability,
            comparator_matrix=comparators,
        )
        summary = build_cross_tissue_specificity_summary(classification)

        iosn = classification[classification["feature"].eq("clr__early_iosn")].iloc[0]
        self.assertEqual(iosn["specificity_class"], "olfactory_specific")
        self.assertEqual(iosn["gateway_age_direction"], "negative")
        self.assertEqual(float(iosn["max_abs_importance"]), 2.0)
        self.assertEqual(int(summary["n_features"].sum()), 3)


if __name__ == "__main__":
    unittest.main()
