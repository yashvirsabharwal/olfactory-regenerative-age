from pathlib import Path
import sys
import tempfile
import unittest

import h5py
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.latent import classify_embedding, inspect_h5ad_obsm, latent_readiness_summary, render_latent_space_plan


class LatentAuditTests(unittest.TestCase):
    def test_classify_embedding_marks_umap_as_visualization_only(self):
        readiness, use = classify_embedding("X_umap", 2)

        self.assertEqual(readiness, "visualization_only")
        self.assertEqual(use, "display_only")

    def test_inspect_h5ad_obsm_reads_shapes_without_matrix(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "toy.h5ad"
            with h5py.File(path, "w") as handle:
                obsm = handle.create_group("obsm")
                obsm.create_dataset("X_umap", data=[[0.0, 1.0], [1.0, 0.0]])
                obsm.create_dataset("X_scvi", data=[[0.0] * 10, [1.0] * 10])

            audit = inspect_h5ad_obsm(path, ["X_scvi", "X_umap"])

        self.assertEqual(audit["embedding_key"].tolist(), ["X_scvi", "X_umap"])
        self.assertEqual(audit.loc[audit["embedding_key"].eq("X_scvi"), "readiness"].iloc[0], "usable_latent")
        self.assertEqual(int(audit.loc[audit["embedding_key"].eq("X_umap"), "n_dimensions"].iloc[0]), 2)

    def test_latent_readiness_summary_requires_non_umap_embedding(self):
        local = pd.DataFrame(
            {
                "embedding_key": ["X_umap"],
                "readiness": ["visualization_only"],
            }
        )
        portal = pd.DataFrame({"portal_embeddings": ["X_umap"]})

        summary = latent_readiness_summary(local, portal)
        plan = render_latent_space_plan(summary, local, portal)

        self.assertEqual(summary.loc[0, "status"], "latent_recompute_required")
        self.assertEqual(summary.loc[0, "usable_local_embeddings"], "none")
        self.assertIn("UMAP alone is not acceptable", plan)


if __name__ == "__main__":
    unittest.main()
