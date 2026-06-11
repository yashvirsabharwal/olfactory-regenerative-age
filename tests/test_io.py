from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))


class AnnDataIoTests(unittest.TestCase):
    def test_inspect_h5ad_uses_backed_mode_on_toy_anndata(self):
        try:
            import anndata as ad
            import numpy as np
            import pandas as pd
            from scipy import sparse
        except ModuleNotFoundError:
            self.skipTest("anndata/scipy are not installed in this runtime")

        from ora.io import inspect_h5ad

        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "toy.h5ad"
            obs = pd.DataFrame(
                {"donor_id": ["d1", "d2"], "sample_id": ["s1", "s2"]},
                index=["c1", "c2"],
            )
            var = pd.DataFrame(index=["g1", "g2", "g3"])
            ad.AnnData(X=sparse.csr_matrix(np.ones((2, 3))), obs=obs, var=var).write_h5ad(path)

            schema, obs_table, var_table = inspect_h5ad(path)

        self.assertEqual(schema["n_obs"], 2)
        self.assertEqual(schema["n_vars"], 3)
        self.assertIn("donor_id", set(obs_table["column"]))
        self.assertEqual(var_table.shape[0], 0)


if __name__ == "__main__":
    unittest.main()

