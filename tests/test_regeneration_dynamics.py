from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from ora.regeneration_dynamics import audit_h5ad_dynamics_inputs, build_dynamics_feasibility


class TestRegenerationDynamics(unittest.TestCase):
    def test_audit_detects_missing_velocity_and_scvi_lineage_support(self):
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "toy.h5ad"
            obs = pd.DataFrame(
                {
                    "fine_celltype": ["Quiescent_HBC", "Early_INP", "Early_iOSN"],
                    "coarse_celltype": ["Resp_HBC", "Olf_INPs", "Olf_iOSNs"],
                    "sample_id": ["s1", "s1", "s2"],
                    "donor_id": ["d1", "d1", "d2"],
                    "flex_version": ["flex_v1", "flex_v1", "flex_v2"],
                    "device_guided": ["T", "T", "F"],
                },
                index=["c1", "c2", "c3"],
            )
            var = pd.DataFrame(index=["g1", "g2"])
            adata = ad.AnnData(X=sparse.csr_matrix(np.ones((3, 2))), obs=obs, var=var)
            adata.obsm["X_scvi"] = np.ones((3, 2))
            adata.write_h5ad(path)

            audit = audit_h5ad_dynamics_inputs([path])
            self.assertFalse(bool(audit.loc[0, "has_spliced"]))
            self.assertFalse(bool(audit.loc[0, "has_unspliced"]))
            self.assertTrue(bool(audit.loc[0, "has_x_scvi"]))
            self.assertTrue(bool(audit.loc[0, "has_lineage_labels"]))

            feasibility = build_dynamics_feasibility(audit)
            status = dict(zip(feasibility["method"], feasibility["status"], strict=True))
            self.assertEqual(status["RNA velocity / scVelo"], "no_go")
            self.assertEqual(status["CellRank velocity kernel"], "no_go")
            self.assertEqual(status["Scanpy diffusion pseudotime"], "feasible_exploratory")


if __name__ == "__main__":
    unittest.main()
