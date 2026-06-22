import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pandas as pd

from ora.neighborhood_programs import score_neighborhood_programs_h5ad


class NeighborhoodProgramTests(unittest.TestCase):
    def test_score_neighborhood_programs_h5ad_scores_memberships(self):
        try:
            import anndata as ad
            from scipy import sparse
        except ModuleNotFoundError:
            self.skipTest("anndata/scipy are not installed in this runtime")

        obs = pd.DataFrame(index=["c0", "c1", "c2", "c3"])
        var = pd.DataFrame({"feature_name": ["TP63", "OMP", "GAP43"]}, index=["g0", "g1", "g2"])
        x = sparse.csr_matrix(np.array([[5, 0, 0], [4, 0, 1], [0, 5, 4], [0, 4, 5]], dtype=float))
        memberships = pd.DataFrame(
            {
                "neighborhood_id": [0, 0, 1, 1],
                "cell_index": [0, 1, 2, 3],
            }
        )
        da = pd.DataFrame(
            {
                "neighborhood_id": [0, 1],
                "top_fine_celltype": ["HBC", "Early_iOSN"],
                "top_coarse_celltype": ["Resp_HBC", "Olf_iOSNs"],
                "age_coef": [-1.0, -0.5],
                "age_fdr": [0.05, 0.2],
                "status": ["tested", "tested"],
            }
        )
        gene_sets = {
            "score": {"log1p": False, "var_symbol_columns": ["feature_name"]},
            "gene_sets": {
                "basal": {"genes": ["TP63"]},
                "immature_neuron": {"genes": ["GAP43", "OMP"]},
            },
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "toy.h5ad"
            ad.AnnData(X=x, obs=obs, var=var).write_h5ad(path)
            scores, summary, coverage = score_neighborhood_programs_h5ad(
                path,
                gene_sets,
                memberships,
                da_table=da,
                run_name="toy",
                chunk_neighborhoods=1,
                log1p=False,
            )

        self.assertEqual(scores.shape[0], 4)
        basal_zero = scores.loc[scores["module"].eq("basal") & scores["neighborhood_id"].eq(0), "program_score"].iloc[0]
        neuron_one = scores.loc[
            scores["module"].eq("immature_neuron") & scores["neighborhood_id"].eq(1),
            "program_score",
        ].iloc[0]
        self.assertAlmostEqual(basal_zero, 4.5)
        self.assertAlmostEqual(neuron_one, 4.5)
        self.assertEqual(int(summary.loc[summary["module"].eq("basal"), "n_significant"].iloc[0]), 1)
        self.assertEqual(set(coverage["module"]), {"basal", "immature_neuron"})


if __name__ == "__main__":
    unittest.main()
