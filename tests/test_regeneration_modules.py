import unittest

import numpy as np
import pandas as pd

from ora.regeneration_modules import (
    build_regeneration_module_age_associations,
    build_regeneration_module_ora_correlations,
    parse_regeneration_module_metadata,
)


class RegenerationModuleTests(unittest.TestCase):
    def test_parse_metadata_preserves_source_fields(self):
        config = {
            "gene_sets": {
                "tp63_hbc_quiescence": {
                    "description": "HBC",
                    "theme": "basal_quiescence",
                    "expected_age_direction": "positive",
                    "source": "literature",
                    "citation": "PMID:1",
                    "genes": ["TP63"],
                }
            }
        }

        metadata = parse_regeneration_module_metadata(config)

        row = metadata.iloc[0]
        self.assertEqual(row["module_feature"], "module_score__tp63_hbc_quiescence")
        self.assertEqual(row["theme"], "basal_quiescence")
        self.assertEqual(row["citation"], "PMID:1")

    def test_age_associations_and_correlations_use_analysis_sets(self):
        donors = [f"d{i}" for i in range(10)]
        module_features = pd.DataFrame(
            {
                "donor_id": donors,
                "module_score__tp63_hbc_quiescence": np.linspace(0, 1, 10),
                "module_score__senescence_sasp": np.linspace(1, 0, 10),
            }
        )
        manifest = pd.DataFrame(
            {
                "donor_id": donors,
                "age": np.linspace(30, 80, 10),
                "usable_for_ora_training": [True] * 10,
                "passes_strict_ora_training_rule": [True] * 8 + [False, False],
            }
        )
        metadata = pd.DataFrame(
            {
                "module": ["tp63_hbc_quiescence"],
                "module_feature": ["module_score__tp63_hbc_quiescence"],
                "theme": ["basal_quiescence"],
                "description": ["HBC"],
                "expected_age_direction": ["positive"],
                "source": ["literature"],
                "citation": ["PMID:1"],
            }
        )
        coverage = pd.DataFrame(
            {
                "module": ["tp63_hbc_quiescence"],
                "n_requested": [1],
                "n_present": [1],
                "coverage_fraction": [1.0],
            }
        )
        ora_features = pd.DataFrame(
            {
                "donor_id": donors,
                "clr__quiescent_hbc": np.linspace(0, 1, 10),
                "prop__goblet": np.linspace(1, 0, 10),
                "module_score__senescence_sasp": np.linspace(0.5, 0.1, 10),
            }
        )
        feature_map = pd.DataFrame(
            {
                "feature": ["clr__quiescent_hbc", "prop__goblet"],
                "primary_theme": ["basal_quiescence", "respiratory_metaplasia_ciliated_goblet"],
                "specificity_class": ["pan_epithelial_regenerative", "airway_nasal_shared"],
            }
        )

        age = build_regeneration_module_age_associations(
            donor_module_features=module_features,
            manifest=manifest,
            module_metadata=metadata,
            coverage=coverage,
        )
        correlations = build_regeneration_module_ora_correlations(
            donor_module_features=module_features,
            ora_feature_matrix=ora_features,
            module_metadata=metadata,
            manifest=manifest,
            feature_map=feature_map,
        )

        primary = age[age["analysis_set"].eq("primary")].iloc[0]
        self.assertEqual(primary["observed_vs_expected"], "aligned")
        self.assertEqual(float(primary["coverage_fraction"]), 1.0)
        self.assertIn("unadjusted", set(age["model_type"]))
        self.assertIn("covariate_adjusted", set(age["model_type"]))
        self.assertIn("strict", set(age["analysis_set"]))
        top = correlations[
            correlations["feature"].eq("clr__quiescent_hbc")
            & correlations["module_feature"].eq("module_score__tp63_hbc_quiescence")
            & correlations["analysis_set"].eq("primary")
        ].iloc[0]
        self.assertAlmostEqual(float(top["pearson_r"]), 1.0)
        self.assertEqual(top["feature_theme"], "basal_quiescence")
        overlap = correlations[
            correlations["module_feature"].eq("module_score__senescence_sasp")
            & correlations["feature"].eq("module_score__senescence_sasp")
            & correlations["analysis_set"].eq("primary")
        ].iloc[0]
        self.assertEqual(overlap["direction"], "positive")


if __name__ == "__main__":
    unittest.main()
