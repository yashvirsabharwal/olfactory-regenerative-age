from pathlib import Path
import importlib.util
import subprocess
import sys
import unittest
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import anndata as ad
import numpy as np
from scipy import sparse

from ora.latent_recompute import (
    latent_recompute_feasibility,
    render_latent_recompute_workflow,
    validate_scvi_pilot,
)


class LatentRecomputeTests(unittest.TestCase):
    def test_latent_recompute_feasibility_records_dependencies_and_scale(self):
        feasibility = latent_recompute_feasibility(
            {"n_obs": 4_000_000, "n_vars": 18_000},
            package_modules={"fakepkg": "module_that_should_not_exist_ora_test"},
            n_top_genes=3000,
            pilot_max_cells=250_000,
        )

        dep = feasibility[feasibility["check"].eq("dependency__fakepkg")].iloc[0]
        cells = feasibility[feasibility["check"].eq("input_cells")].iloc[0]
        self.assertEqual(dep["status"], "missing")
        self.assertEqual(cells["status"], "large")
        self.assertIn("pilot subset", cells["recommendation"])

    def test_render_latent_recompute_workflow_mentions_missing_dependencies(self):
        feasibility = latent_recompute_feasibility(
            {"n_obs": 1000, "n_vars": 5000},
            package_modules={"fakepkg": "module_that_should_not_exist_ora_test"},
        )

        text = render_latent_recompute_workflow(
            feasibility,
            h5ad_path="data/raw/gateway.h5ad",
            output_h5ad="data/processed/gateway_scvi_pilot.h5ad",
        )

        self.assertIn("Missing latent dependencies: fakepkg", text)
        self.assertIn("scripts/run_scvi_latent.py", text)

    def test_scvi_runner_reads_bounded_backed_subset(self):
        script_path = Path(__file__).resolve().parents[1] / "scripts" / "run_scvi_latent.py"
        spec = importlib.util.spec_from_file_location("run_scvi_latent_for_test", script_path)
        module = importlib.util.module_from_spec(spec)
        assert spec.loader is not None
        spec.loader.exec_module(module)
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "toy.h5ad"
            adata = ad.AnnData(X=sparse.csr_matrix(np.arange(60).reshape(10, 6)))
            adata.write_h5ad(path)

            subset = module._read_pilot_h5ad(ad, str(path), max_cells=4, seed=1)

        self.assertEqual(subset.n_obs, 4)
        self.assertEqual(subset.n_vars, 6)

    def test_scvi_runner_reads_all_cells_with_selected_genes(self):
        script_path = Path(__file__).resolve().parents[1] / "scripts" / "run_scvi_latent.py"
        spec = importlib.util.spec_from_file_location("run_scvi_latent_for_test", script_path)
        module = importlib.util.module_from_spec(spec)
        assert spec.loader is not None
        spec.loader.exec_module(module)
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "toy.h5ad"
            adata = ad.AnnData(X=np.arange(60).reshape(10, 6))
            adata.var_names = [f"gene_{idx}" for idx in range(6)]
            adata.write_h5ad(path)

            subset = module._read_pilot_h5ad(
                ad,
                str(path),
                max_cells=None,
                seed=1,
                selected_var_names=("gene_1", "gene_3"),
            )

        self.assertEqual(subset.n_obs, 10)
        self.assertEqual(subset.var_names.tolist(), ["gene_1", "gene_3"])

    def test_scvi_runner_can_stratify_backed_subset(self):
        script_path = Path(__file__).resolve().parents[1] / "scripts" / "run_scvi_latent.py"
        spec = importlib.util.spec_from_file_location("run_scvi_latent_for_test", script_path)
        module = importlib.util.module_from_spec(spec)
        assert spec.loader is not None
        spec.loader.exec_module(module)
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "toy.h5ad"
            adata = ad.AnnData(X=np.arange(120).reshape(20, 6))
            adata.obs["group"] = ["rare"] * 2 + ["common"] * 18
            adata.write_h5ad(path)

            subset = module._read_pilot_h5ad(
                ad,
                str(path),
                max_cells=4,
                seed=1,
                sampling_strategy="stratified",
                stratify_keys=("group",),
            )

        self.assertEqual(subset.obs["group"].tolist().count("rare"), 2)

    def test_build_reduced_h5ad_writes_selected_genes_in_chunks(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            source = tmp / "source.h5ad"
            genes = tmp / "genes.txt"
            out = tmp / "reduced.h5ad"
            chunks = tmp / "chunks"
            adata = ad.AnnData(X=np.arange(60).reshape(10, 6))
            adata.var_names = [f"gene_{idx}" for idx in range(6)]
            adata.obs["sample_id"] = [f"sample_{idx % 2}" for idx in range(10)]
            adata.write_h5ad(source)
            genes.write_text("gene_1\ngene_3\ngene_5\n", encoding="utf-8")

            subprocess.run(
                [
                    sys.executable,
                    "scripts/build_reduced_h5ad.py",
                    "--h5ad",
                    str(source),
                    "--gene-list-file",
                    str(genes),
                    "--out",
                    str(out),
                    "--chunk-dir",
                    str(chunks),
                    "--chunk-size",
                    "4",
                ],
                check=True,
                cwd=Path(__file__).resolve().parents[1],
            )

            reduced = ad.read_h5ad(out)

        self.assertEqual(reduced.n_obs, 10)
        self.assertEqual(reduced.var_names.tolist(), ["gene_1", "gene_3", "gene_5"])
        self.assertEqual(reduced.obs["sample_id"].tolist()[:3], ["sample_0", "sample_1", "sample_0"])

    def test_validate_scvi_pilot_reports_embedding_and_sparse_metadata(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "pilot.h5ad"
            adata = ad.AnnData(X=np.arange(60).reshape(10, 6))
            adata.obsm["X_scvi"] = np.ones((10, 10))
            adata.obs["sample_id"] = ["rare"] + ["common"] * 9
            adata.obs["donor_id"] = ["d1"] * 5 + ["d2"] * 5
            adata.write_h5ad(path)

            validation = validate_scvi_pilot(
                path,
                metadata_columns=("sample_id", "donor_id", "missing_column"),
            )

        dim = validation[validation["check"].eq("embedding_dimensions")].iloc[0]
        sample = validation[validation["check"].eq("metadata__sample_id")].iloc[0]
        missing = validation[validation["check"].eq("metadata__missing_column")].iloc[0]
        self.assertEqual(dim["status"], "ok")
        self.assertEqual(sample["status"], "sparse_levels")
        self.assertEqual(missing["status"], "missing")


if __name__ == "__main__":
    unittest.main()
