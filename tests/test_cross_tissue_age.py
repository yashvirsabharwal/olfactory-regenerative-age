import unittest

import numpy as np
import pandas as pd

from ora.cross_tissue_age import (
    CONTEXT_SCOPE,
    PRIMARY_SCOPE,
    estimate_cross_tissue_age_effects,
    harmonize_cell_group,
    map_ora_feature_to_cross_tissue_feature,
    parse_cellxgene_age,
)


class CrossTissueAgeEffectTests(unittest.TestCase):
    def test_parse_cellxgene_age_handles_years_and_developmental_weeks(self):
        self.assertEqual(parse_cellxgene_age("70-year-old stage"), 70.0)
        self.assertAlmostEqual(parse_cellxgene_age("31st week post-fertilization stage"), 31 / 52)
        self.assertTrue(np.isnan(parse_cellxgene_age("unknown")))

    def test_harmonize_cell_group_maps_airway_and_immune_labels(self):
        self.assertEqual(harmonize_cell_group("nasal mucosa goblet cell"), "goblet")
        self.assertEqual(harmonize_cell_group("respiratory tract multiciliated cell"), "ciliated")
        self.assertEqual(harmonize_cell_group("Neutrophils & macrophages"), "myeloid")
        self.assertEqual(harmonize_cell_group("professional antigen presenting cell"), "antigen_presenting")

    def test_map_ora_feature_to_cross_tissue_feature(self):
        self.assertEqual(map_ora_feature_to_cross_tissue_feature("prop__goblet"), ["cell_group__goblet"])
        self.assertEqual(map_ora_feature_to_cross_tissue_feature("clr__cdc1"), ["cell_group__antigen_presenting"])
        self.assertEqual(
            map_ora_feature_to_cross_tissue_feature("module_score__senescence_sasp"),
            ["module_score__senescence_sasp"],
        )
        self.assertEqual(map_ora_feature_to_cross_tissue_feature("ratio__mature_mosn_to_iosn"), [])

    def test_estimate_age_effects_uses_adult_primary_scope(self):
        donor_features = pd.DataFrame(
            {
                "dataset_id": ["toy"] * 6,
                "dataset_title": ["Toy"] * 6,
                "tissue_class": ["nasal"] * 6,
                "donor_id": [f"d{i}" for i in range(6)],
                "age_years": [8, 12, 30, 40, 50, 60],
                "age_context": ["child_or_developmental", "child_or_developmental", "adult", "adult", "adult", "adult"],
                "n_cells": [100] * 6,
                "cell_group__goblet": [0.01, 0.02, 0.1, 0.2, 0.3, 0.4],
                "module_score__senescence_sasp": [0.2, 0.2, 1.0, 1.1, 1.2, 1.3],
            }
        )

        effects = estimate_cross_tissue_age_effects(donor_features, min_donors=4)

        adult_goblet = effects[
            effects["analysis_scope"].eq(PRIMARY_SCOPE) & effects["feature"].eq("cell_group__goblet")
        ].iloc[0]
        context_goblet = effects[
            effects["analysis_scope"].eq(CONTEXT_SCOPE) & effects["feature"].eq("cell_group__goblet")
        ].iloc[0]
        self.assertEqual(adult_goblet["status"], "ok")
        self.assertEqual(adult_goblet["direction"], "positive")
        self.assertEqual(context_goblet["status"], "context_only")


if __name__ == "__main__":
    unittest.main()
