from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.modules import parse_gene_sets, resolve_gene_sets, score_gene_sets_h5ad


def gateway_like_config():
    return {
        "columns": {
            "donor_id": ["donor_id"],
            "sample_id": ["sample_id"],
            "age": ["age"],
            "disease": ["disease"],
            "coarse_cell_type": ["coarse_cell_type"],
            "fine_cell_type": ["fine_cell_type"],
            "sex": ["sex"],
        },
        "healthy_values": ["healthy"],
        "ndd_values": {"ad": ["ad"], "pd": ["pd"]},
    }


class ModuleTests(unittest.TestCase):
    def test_parse_gene_sets_accepts_dict_entries(self):
        parsed = parse_gene_sets({"gene_sets": {"m": {"description": "module", "genes": ["A", "B"]}}})

        self.assertEqual(parsed[0].name, "m")
        self.assertEqual(parsed[0].genes, ("A", "B"))
        self.assertEqual(parsed[0].description, "module")

    def test_resolve_gene_sets_uses_feature_name(self):
        var = pd.DataFrame({"feature_name": ["TP63", "OMP", "SNCA"]}, index=["ENSG1", "ENSG2", "ENSG3"])
        gene_sets = parse_gene_sets({"gene_sets": {"olf": {"genes": ["TP63", "OMP", "MISSING"]}}})

        resolved, coverage = resolve_gene_sets(var, var.index, gene_sets)

        self.assertEqual(resolved["olf"], [0, 1])
        self.assertEqual(int(coverage.loc[0, "n_present"]), 2)
        self.assertEqual(coverage.loc[0, "missing_genes"], "MISSING")

    def test_score_gene_sets_h5ad_chunks_and_aggregates(self):
        try:
            import anndata as ad
            from scipy import sparse
        except ModuleNotFoundError:
            self.skipTest("anndata/scipy are not installed in this runtime")

        obs = pd.DataFrame(
            {
                "donor_id": ["d1", "d1", "d2"],
                "sample_id": ["s1", "s1", "s2"],
                "age": [40, 40, 70],
                "disease": ["healthy", "healthy", "healthy"],
                "sex": ["female", "female", "male"],
                "coarse_cell_type": ["HBC", "Neuron", "Neuron"],
                "fine_cell_type": ["qHBC", "mOSN", "mOSN"],
            },
            index=["c1", "c2", "c3"],
        )
        var = pd.DataFrame({"feature_name": ["TP63", "OMP", "SNCA"]}, index=["ENSG1", "ENSG2", "ENSG3"])
        x = sparse.csr_matrix(np.array([[0, 1, 4], [2, 3, 0], [4, 5, 6]], dtype=float))
        gene_sets = {"gene_sets": {"basal": {"genes": ["TP63"]}, "mixed": {"genes": ["TP63", "SNCA"]}}}

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "toy.h5ad"
            ad.AnnData(X=x, obs=obs, var=var).write_h5ad(path)

            result = score_gene_sets_h5ad(
                path,
                gateway_like_config(),
                gene_sets,
                chunk_size=2,
                log1p=False,
            )

        self.assertEqual(set(result.coverage["module"]), {"basal", "mixed"})
        self.assertIn("module_score__mixed", set(result.donor_features.columns))
        d1_mixed = result.donor_features.loc[result.donor_features["donor_id"].eq("d1"), "module_score__mixed"].iloc[0]
        self.assertAlmostEqual(d1_mixed, 1.5)
        self.assertEqual(int(result.summary["n_cells"].sum()), 6)


if __name__ == "__main__":
    unittest.main()
