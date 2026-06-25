from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.feature_ablation import (
    build_feature_family_matrices,
    build_technical_covariate_matrix,
    maybe_build_pseudobulk_pc_features,
    run_feature_family_ablation,
)


class FeatureFamilyAblationTests(unittest.TestCase):
    def test_build_feature_family_matrices_splits_expected_families(self):
        features = pd.DataFrame(
            {
                "donor_id": ["d1", "d2"],
                "prop__hbc": [0.1, 0.2],
                "clr__hbc": [0.0, 1.0],
                "ratio__inp_to_hbc": [1.0, 2.0],
                "module_score__stress": [0.5, 0.4],
                "scvi_global_mean__dim01": [1.0, 2.0],
                "scvi_state_mean__hbc__dim01": [3.0, 4.0],
            }
        )
        manifest = pd.DataFrame({"donor_id": ["d1", "d2"], "sample_id": ["s1", "s2"], "sex": ["F", "M"]})

        matrices, feasibility = build_feature_family_matrices(features, manifest)

        self.assertEqual(matrices["proportions_only"].columns.tolist(), ["donor_id", "prop__hbc"])
        self.assertEqual(matrices["modules_only"].columns.tolist(), ["donor_id", "module_score__stress"])
        self.assertIn("ora_scvi_hybrid", matrices)
        self.assertIn("technical_covariates_only", matrices)
        self.assertEqual(feasibility.loc[feasibility["feature_set"].eq("pseudobulk_expression_pcs"), "status"].iloc[0], "skipped_missing_inputs")

    def test_build_technical_covariate_matrix_one_hot_encodes_manifest(self):
        manifest = pd.DataFrame(
            {
                "donor_id": ["d1", "d2"],
                "sample_id": ["s1", "s2"],
                "sex": ["F", "M"],
                "chemistry": ["v1", "v2"],
                "total_cells": [99, 999],
            }
        )

        technical = build_technical_covariate_matrix(manifest)

        self.assertIn("technical__log10_total_cells", technical.columns)
        self.assertIn("technical__sex_F", technical.columns)
        self.assertIn("technical__chemistry_v2", technical.columns)

    def test_maybe_build_pseudobulk_pc_features_aggregates_to_donors(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            counts = tmp / "counts.tsv"
            metadata = tmp / "metadata.tsv"
            pd.DataFrame(
                {
                    "gene_id": [f"g{i}" for i in range(5)],
                    "PB1": [10, 0, 5, 1, 0],
                    "PB2": [0, 10, 5, 0, 1],
                    "PB3": [5, 5, 0, 10, 0],
                }
            ).to_csv(counts, sep="\t", index=False)
            pd.DataFrame(
                {
                    "pseudobulk_id": ["PB1", "PB2", "PB3"],
                    "donor_id": ["d1", "d2", "d3"],
                }
            ).to_csv(metadata, sep="\t", index=False)

            pcs = maybe_build_pseudobulk_pc_features(counts, metadata, n_components=2, n_variable_genes=4)

        assert pcs is not None
        self.assertEqual(pcs.shape, (3, 3))
        self.assertEqual(pcs.columns.tolist(), ["donor_id", "pseudobulk_pc__PC01", "pseudobulk_pc__PC02"])

    def test_run_feature_family_ablation_outputs_summary_and_deltas(self):
        donors = [f"d{i}" for i in range(8)]
        manifest = pd.DataFrame(
            {
                "donor_id": donors,
                "sample_id": donors,
                "age": np.linspace(30, 70, 8),
                "sex": ["F", "M"] * 4,
                "chemistry": ["v2"] * 8,
                "collection_method": ["device"] * 8,
                "site": ["site"] * 8,
                "disease_group": ["healthy"] * 8,
                "usable_for_ora_training": [True] * 8,
                "total_cells": np.arange(8) + 100,
            }
        )
        features = pd.DataFrame(
            {
                "donor_id": donors,
                "prop__hbc": np.linspace(0, 1, 8),
                "clr__hbc": np.linspace(1, 0, 8),
                "ratio__inp_to_hbc": np.linspace(0.2, 0.8, 8),
                "module_score__stress": np.linspace(0.1, 0.9, 8),
                "scvi_global_mean__dim01": np.linspace(0.3, 0.7, 8),
                "scvi_state_mean__hbc__dim01": np.linspace(0.5, 0.1, 8),
            }
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            result = run_feature_family_ablation(
                feature_matrix=features,
                manifest=manifest,
                model_config={"outer_cv_folds": 2, "random_seed": 1, "model_names": ["ridge"]},
                output_dir=tmpdir,
                models=["ridge"],
                repeats=1,
                n_permutations=0,
            )

        self.assertFalse(result.summary.empty)
        self.assertFalse(result.deltas.empty)
        self.assertIn("calibration_slope_mean", result.summary.columns)


if __name__ == "__main__":
    unittest.main()
