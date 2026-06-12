from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.pseudobulk import aggregate_targeted_pseudobulk_h5ad, run_covariate_pseudobulk_de, run_pseudobulk_de


def gateway_like_config():
    return {
        "columns": {
            "donor_id": ["donor_id"],
            "sample_id": ["sample_id"],
            "age": ["age"],
            "disease": ["disease"],
            "coarse_cell_type": ["coarse_cell_type"],
            "fine_cell_type": ["fine_cell_type"],
            "n_counts": ["nCount_RNA"],
        },
        "healthy_values": ["healthy"],
        "ndd_values": {"ad": ["ad"], "pd": ["pd"]},
    }


class PseudobulkTests(unittest.TestCase):
    def test_targeted_pseudobulk_aggregates_feature_name_genes(self):
        try:
            import anndata as ad
            from scipy import sparse
        except ModuleNotFoundError:
            self.skipTest("anndata/scipy are not installed in this runtime")

        obs = pd.DataFrame(
            {
                "donor_id": ["h1", "h1", "a1"],
                "sample_id": ["s1", "s1", "s2"],
                "age": [40, 40, 75],
                "disease": ["healthy", "healthy", "ad"],
                "coarse_cell_type": ["HBC", "HBC", "HBC"],
                "fine_cell_type": ["qHBC", "qHBC", "qHBC"],
                "nCount_RNA": [10, 20, 30],
            },
            index=["c1", "c2", "c3"],
        )
        var = pd.DataFrame({"feature_name": ["TP63", "OMP", "SNCA"]}, index=["ENSG1", "ENSG2", "ENSG3"])
        x = sparse.csr_matrix(np.array([[1, 0, 3], [2, 5, 0], [4, 6, 8]], dtype=float))

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "toy.h5ad"
            ad.AnnData(X=x, obs=obs, var=var).write_h5ad(path)
            result = aggregate_targeted_pseudobulk_h5ad(
                path,
                gateway_like_config(),
                ["TP63", "SNCA", "MISSING"],
                chunk_size=2,
                min_donors=1,
            )

        self.assertEqual(int(result.coverage.loc[0, "n_present"]), 2)
        self.assertEqual(result.metadata["n_cells"].sum(), 3)
        h1_tp63 = result.counts[
            result.counts["donor_id"].eq("h1") & result.counts["gene"].eq("TP63")
        ]["count"].iloc[0]
        self.assertEqual(h1_tp63, 3)

    def test_pseudobulk_de_reports_ok_contrast(self):
        counts = np.array(
            [
                [100, 10],
                [110, 10],
                [120, 10],
                [10, 80],
                [12, 90],
                [11, 85],
            ],
            dtype=float,
        )
        group_meta = pd.DataFrame(
            {
                "donor_id": ["a1", "a2", "a3", "h1", "h2", "h3"],
                "disease_group": ["ad", "ad", "ad", "healthy", "healthy", "healthy"],
                "fine_cell_type": ["qHBC"] * 6,
                "n_cells": [10] * 6,
                "sum_n_counts": counts.sum(axis=1),
            }
        )

        de = run_pseudobulk_de(
            counts,
            group_meta,
            ["TP63", "OMP"],
            contrasts=[("ad", "healthy")],
            min_donors=3,
        )

        tp63 = de[de["gene"].eq("TP63")].iloc[0]
        self.assertEqual(tp63["status"], "ok")
        self.assertGreater(tp63["log2fc"], 0)

    def test_covariate_pseudobulk_de_reports_adjusted_effect(self):
        metadata = pd.DataFrame(
            {
                "donor_id": ["a1", "a2", "a3", "h1", "h2", "h3"],
                "sample_id": ["s1", "s2", "s3", "s4", "s5", "s6"],
                "disease_group": ["ad", "ad", "ad", "healthy", "healthy", "healthy"],
                "coarse_cell_type": ["HBC"] * 6,
                "fine_cell_type": ["qHBC"] * 6,
                "n_cells": [10] * 6,
                "sum_n_counts": [1000] * 6,
            }
        )
        counts = pd.DataFrame(
            {
                "donor_id": ["a1", "a2", "a3", "h1", "h2", "h3"],
                "sample_id": ["s1", "s2", "s3", "s4", "s5", "s6"],
                "disease_group": ["ad", "ad", "ad", "healthy", "healthy", "healthy"],
                "coarse_cell_type": ["HBC"] * 6,
                "fine_cell_type": ["qHBC"] * 6,
                "gene": ["TP63"] * 6,
                "count": [90, 100, 110, 10, 12, 14],
            }
        )
        manifest = pd.DataFrame(
            {
                "donor_id": ["a1", "a2", "a3", "h1", "h2", "h3"],
                "age": [40, 50, 60, 40, 50, 60],
                "sex": ["F", "M", "F", "F", "M", "F"],
                "chemistry": ["v2"] * 6,
                "collection_method": ["device"] * 6,
            }
        )

        de = run_covariate_pseudobulk_de(
            counts,
            metadata,
            manifest,
            genes=["TP63"],
            contrasts=[("ad", "healthy")],
            covariates=["age", "sex", "chemistry", "collection_method"],
            min_donors=3,
        )

        row = de.iloc[0]
        self.assertEqual(row["status"], "ok")
        self.assertTrue(np.isfinite(row["p_value"]))
        self.assertGreater(row["log2fc_adjusted"], 0)
        self.assertIn("age", row["covariates"])


if __name__ == "__main__":
    unittest.main()
