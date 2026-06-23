from pathlib import Path
import gzip
import sys
import tarfile
import tempfile
import unittest

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.external import (
    build_external_10x_marker_mapped_anndata,
    external_candidate_matrix,
    external_dataset_summary,
    external_mapped_feature_concordance,
    external_marker_age_concordance,
    external_module_contrasts,
    external_validation_evidence_summary,
    feature_matrix_contract_summary,
    inspect_external_archive,
    parse_published_gene_lists,
    parse_geo_series_matrix_metadata,
    published_gene_list_coverage,
    score_external_10x_marker_composition,
    score_external_10x_modules,
    validate_external_feature_matrix,
)


class ExternalValidationTests(unittest.TestCase):
    def test_external_dataset_summary_tracks_file_readiness(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            feature_matrix = root / "features.tsv"
            feature_matrix.write_text("donor_id\tage\nD1\t40\n", encoding="utf-8")
            config = {
                "datasets": {
                    "toy": {
                        "title": "Toy validation",
                        "status": "configured",
                        "validation_use": "unit test",
                        "species": "human",
                        "tissue": "olfactory",
                        "disease_context": ["healthy"],
                        "expected_level": "donor_feature_matrix",
                        "required_files": {"feature_matrix": "features.tsv", "expression": None, "metadata": None},
                    }
                }
            }

            summary = external_dataset_summary(config, base_dir=root)

        self.assertEqual(summary.loc[0, "dataset_id"], "toy")
        self.assertTrue(bool(summary.loc[0, "ready_for_feature_validation"]))
        self.assertFalse(bool(summary.loc[0, "ready_for_raw_adapter"]))
        self.assertIn("metadata", summary.loc[0, "files_missing"])
        self.assertEqual(summary.loc[0, "readiness_class"], "ready_feature_matrix")

    def test_external_candidate_matrix_classifies_direct_and_context_sources(self):
        config = {
            "datasets": {
                "direct": {
                    "title": "Olfactory aging",
                    "accession": "GSE184117",
                    "status": "download_ready",
                    "validation_use": "aging feature replication",
                    "species": "human",
                    "tissue": "olfactory epithelium",
                    "disease_context": ["aging", "presbyosmia"],
                    "expected_level": "single_cell",
                    "required_files": {"expression": None, "metadata": None, "feature_matrix": None},
                },
                "bulk": {
                    "title": "Bulk context",
                    "status": "download_ready_bulk_reference",
                    "validation_use": "marker sanity",
                    "species": "human",
                    "tissue": "olfactory neuroepithelium",
                    "disease_context": ["healthy"],
                    "expected_level": "bulk_rnaseq",
                    "required_files": {"expression": None, "metadata": None, "feature_matrix": None},
                },
                "context": {
                    "title": "Long COVID olfactory",
                    "status": "download_ready_context",
                    "validation_use": "disease context",
                    "species": "human",
                    "tissue": "olfactory epithelium",
                    "disease_context": ["long_covid_hyposmia", "normosmic_control"],
                    "expected_level": "single_cell",
                    "required_files": {"expression": None, "metadata": None, "feature_matrix": None},
                },
            }
        }

        matrix = external_candidate_matrix(config)

        direct = matrix[matrix["dataset_id"].eq("direct")].iloc[0]
        bulk = matrix[matrix["dataset_id"].eq("bulk")].iloc[0]
        self.assertEqual(direct["validation_class"], "direct_small_n_mapping_candidate")
        self.assertEqual(direct["supports_primary_claim"], "guarded_small_n_only")
        self.assertEqual(bulk["validation_class"], "context_only_bulk_marker")
        self.assertEqual(bulk["supports_primary_claim"], "no")
        context = matrix[matrix["dataset_id"].eq("context")].iloc[0]
        self.assertEqual(context["validation_class"], "human_olfactory_context_adapter_candidate")
        self.assertEqual(context["adapter_status"], "download_available_missing_local_files")

    def test_inspect_external_archive_classifies_10x_members(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            tar_path = root / "raw.tar"
            for name in ["S1_matrix.mtx", "S1_barcodes.tsv", "S1_features.tsv"]:
                (root / name).write_text("x\n", encoding="utf-8")
            with tarfile.open(tar_path, "w") as archive:
                for name in ["S1_matrix.mtx", "S1_barcodes.tsv", "S1_features.tsv"]:
                    archive.add(root / name, arcname=name)

            inventory = inspect_external_archive(tar_path, dataset_id="toy")

        self.assertEqual(set(inventory["role"]), {"matrix", "barcodes", "features"})
        self.assertEqual(set(inventory["sample_guess"]), {"S1"})

    def test_validate_external_feature_matrix_reports_harmonization(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            feature_matrix = root / "features.tsv"
            feature_matrix.write_text(
                "donor_id\tage\tdisease_group\tprop__hbc\tmodule_score__regen\n"
                "d1\t40\thealthy\t0.1\t1.2\n",
                encoding="utf-8",
            )
            gateway = pd.DataFrame(columns=["donor_id", "prop__hbc", "prop__mosn", "module_score__regen"])
            config = {
                "feature_matrix_contract": {
                    "required_columns": ["donor_id", "age", "disease_group"],
                    "optional_covariates": ["sex"],
                    "accepted_feature_prefixes": ["prop__", "module_score__"],
                }
            }

            summary, harmonization = validate_external_feature_matrix(
                feature_matrix,
                config,
                gateway_features=gateway,
                dataset_id="toy",
            )

        self.assertEqual(summary.loc[0, "status"], "ok")
        self.assertEqual(int(summary.loc[0, "feature_columns"]), 2)
        missing = harmonization[harmonization["status"].eq("missing_from_external")]
        self.assertEqual(missing["feature"].tolist(), ["prop__mosn"])

    def test_published_gene_list_coverage_uses_feature_name(self):
        config = {
            "published_gene_lists": {
                "regen": {"description": "Regeneration", "genes": ["TP63", "OMP", "MISSING"]}
            }
        }
        var = pd.DataFrame({"feature_name": ["TP63", "OMP"]}, index=["ENSG1", "ENSG2"])

        coverage = published_gene_list_coverage(config, var, var.index, ["feature_name"])

        self.assertEqual(coverage.loc[0, "gene_list"], "regen")
        self.assertEqual(int(coverage.loc[0, "n_present"]), 2)
        self.assertEqual(coverage.loc[0, "missing_genes"], "MISSING")

    def test_parse_gene_lists_and_contract_summary(self):
        config = {
            "published_gene_lists": {"aging": {"genes": ["CDKN1A", "CDKN2A"]}},
            "feature_matrix_contract": {
                "required_columns": ["donor_id", "age"],
                "optional_covariates": ["sex"],
                "accepted_feature_prefixes": ["prop__"],
            },
        }

        gene_lists = parse_published_gene_lists(config)
        contract = feature_matrix_contract_summary(config)

        self.assertEqual(gene_lists[0].genes, ("CDKN1A", "CDKN2A"))
        self.assertEqual(set(contract["kind"]), {"required_column", "optional_covariate", "accepted_feature_prefix"})

    def test_parse_geo_series_matrix_metadata_normalizes_samples(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "series.txt.gz"
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write(
                    '!Sample_title\t"Control"\t"Presbyosmic"\t"Culture"\n'
                    '!Sample_geo_accession\t"GSM1"\t"GSM2"\t"GSM3"\n'
                    '!Sample_characteristics_ch1\t"age: 71"\t"age: 78"\t"age: 51"\n'
                    '!Sample_characteristics_ch1\t"disease state: Normosmia"\t'
                    '"disease state: Moderate Hyposmia"\t"disease state: Normosmia"\n'
                    '!Sample_description\t"Control 1"\t"Presbyosmic 1"\t"Culture"\n'
                    '!Sample_supplementary_file_1\t'
                    '"ftp://x/GSM1_H2020_1_Control.barcodes.tsv.gz"\t'
                    '"ftp://x/GSM2_H2020_2_Presbyosmic.barcodes.tsv.gz"\t'
                    '"ftp://x/GSM3_H2020_3_culture.barcodes.tsv.gz"\n'
                    "!series_matrix_table_begin\n"
                )

            metadata = parse_geo_series_matrix_metadata(path, dataset_id="toy")

        self.assertEqual(metadata["disease_group"].tolist(), ["healthy", "presbyosmia", "culture"])
        self.assertEqual(metadata["donor_id"].tolist(), ["H2020_1", "H2020_2", "H2020_3"])
        self.assertEqual(metadata["usable_for_external_validation"].tolist(), [True, True, False])

    def test_score_external_10x_modules_and_contrasts(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            metadata = pd.DataFrame(
                {
                    "dataset_id": ["toy", "toy"],
                    "sample_id": ["GSM1", "GSM2"],
                    "donor_id": ["d1", "d2"],
                    "sample_prefix": ["GSM1_Control", "GSM2_Presbyosmic"],
                    "age": [70, 78],
                    "disease_state": ["Normosmia", "Moderate Hyposmia"],
                    "disease_group": ["healthy", "presbyosmia"],
                    "sample_class": ["biopsy", "biopsy"],
                }
            )
            for prefix, counts in {
                "GSM1_Control": [(1, 1, 5), (2, 1, 1), (2, 2, 1)],
                "GSM2_Presbyosmic": [(1, 1, 2), (2, 1, 5), (2, 2, 4)],
            }.items():
                _write_gzip(root / f"{prefix}.features.tsv.gz", "ENSG1\tTP63\tGene Expression\nENSG2\tKRT5\tGene Expression\n")
                _write_gzip(root / f"{prefix}.barcodes.tsv.gz", "c1\nc2\n")
                matrix_lines = ["%%MatrixMarket matrix coordinate integer general\n", "%\n", "2 2 3\n"]
                matrix_lines.extend(f"{i} {j} {x}\n" for i, j, x in counts)
                _write_gzip(root / f"{prefix}.matrix.mtx.gz", "".join(matrix_lines))
            tar_path = root / "toy.tar"
            with tarfile.open(tar_path, "w") as archive:
                for path in root.glob("*.gz"):
                    archive.add(path, arcname=path.name)

            qc, scores, contrasts = score_external_10x_modules(
                tar_path,
                metadata,
                {"gene_sets": {"basal": {"genes": ["TP63", "KRT5"]}}},
                dataset_id="toy",
            )

        self.assertEqual(qc.shape[0], 2)
        self.assertEqual(scores["module"].unique().tolist(), ["basal"])
        self.assertEqual(contrasts.loc[0, "module"], "basal")
        self.assertEqual(int(contrasts.loc[0, "n_healthy"]), 1)

    def test_score_external_10x_marker_composition_assigns_marker_panels(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            metadata = pd.DataFrame(
                {
                    "dataset_id": ["toy", "toy"],
                    "sample_id": ["GSM1", "GSM2"],
                    "donor_id": ["d1", "d2"],
                    "sample_prefix": ["GSM1_Control", "GSM2_Presbyosmic"],
                    "age": [70, 78],
                    "disease_state": ["Normosmia", "Moderate Hyposmia"],
                    "disease_group": ["healthy", "presbyosmia"],
                    "sample_class": ["biopsy", "biopsy"],
                }
            )
            for prefix, counts in {
                "GSM1_Control": [(1, 1, 10), (1, 2, 8), (2, 3, 9)],
                "GSM2_Presbyosmic": [(1, 1, 9), (2, 2, 12), (2, 3, 11)],
            }.items():
                _write_gzip(
                    root / f"{prefix}.features.tsv.gz",
                    "ENSG1\tTP63\tGene Expression\nENSG2\tOMP\tGene Expression\n",
                )
                _write_gzip(root / f"{prefix}.barcodes.tsv.gz", "c1\nc2\nc3\n")
                matrix_lines = ["%%MatrixMarket matrix coordinate integer general\n", "%\n", "2 3 3\n"]
                matrix_lines.extend(f"{i} {j} {x}\n" for i, j, x in counts)
                _write_gzip(root / f"{prefix}.matrix.mtx.gz", "".join(matrix_lines))
            tar_path = root / "toy.tar"
            with tarfile.open(tar_path, "w") as archive:
                for path in root.glob("*.gz"):
                    archive.add(path, arcname=path.name)

            coverage, composition, contrasts = score_external_10x_marker_composition(
                tar_path,
                metadata,
                {"hbc": {"genes": ["TP63"]}, "mature_osn": {"genes": ["OMP"]}},
                dataset_id="toy",
            )

        self.assertEqual(set(coverage["marker_panel"]), {"hbc", "mature_osn"})
        hbc_control = composition[
            composition["sample_id"].eq("GSM1") & composition["marker_panel"].eq("hbc")
        ].iloc[0]
        self.assertEqual(int(hbc_control["n_cells_assigned"]), 2)
        self.assertIn("mature_osn", contrasts["marker_panel"].tolist())

    def test_build_external_10x_marker_mapped_anndata_writes_feature_contract(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            metadata = pd.DataFrame(
                {
                    "dataset_id": ["toy"],
                    "sample_id": ["GSM1"],
                    "donor_id": ["d1"],
                    "sample_prefix": ["GSM1_Control"],
                    "age": [70],
                    "disease_state": ["Normosmia"],
                    "disease_group": ["healthy"],
                    "sample_class": ["biopsy"],
                    "usable_for_external_validation": [True],
                    "chemistry": ["10x"],
                    "collection_method": ["biopsy"],
                }
            )
            _write_gzip(
                root / "GSM1_Control.features.tsv.gz",
                "ENSG1\tTP63\tGene Expression\nENSG2\tOMP\tGene Expression\n",
            )
            _write_gzip(root / "GSM1_Control.barcodes.tsv.gz", "c1\nc2\nc3\n")
            matrix_lines = ["%%MatrixMarket matrix coordinate integer general\n", "%\n", "2 3 3\n"]
            matrix_lines.extend(["1 1 10\n", "2 2 12\n", "2 3 11\n"])
            _write_gzip(root / "GSM1_Control.matrix.mtx.gz", "".join(matrix_lines))
            tar_path = root / "toy.tar"
            with tarfile.open(tar_path, "w") as archive:
                for path in root.glob("*.gz"):
                    archive.add(path, arcname=path.name)

            adata, qc, donor_features = build_external_10x_marker_mapped_anndata(
                tar_path,
                metadata,
                {"quiescent_hbc": {"genes": ["TP63"]}, "mature_osn": {"genes": ["OMP"]}},
                dataset_id="toy",
            )

        self.assertEqual(adata.n_obs, 3)
        self.assertEqual(set(adata.obs["mapped_cell_state"]), {"quiescent_hbc", "mature_osn"})
        self.assertEqual(qc.loc[0, "status"], "marker_reference_mapped")
        self.assertIn("prop__mature_osn", donor_features.columns)
        self.assertIn("clr__quiescent_hbc", donor_features.columns)
        self.assertIn("ratio__mature_osn_to_quiescent_hbc", donor_features.columns)

    def test_external_mapped_feature_concordance_compares_gateway_age_direction(self):
        mapped = pd.DataFrame(
            {
                "disease_group": ["healthy", "healthy", "presbyosmia", "presbyosmia"],
                "prop__mature_osn": [0.1, 0.2, 0.4, 0.5],
            }
        )
        age = pd.DataFrame(
            {
                "feature": ["prop__fully_mature_mosn"],
                "direction": ["positive"],
                "beta_per_10_years": [0.2],
                "p_value": [0.01],
                "fdr": [0.05],
            }
        )

        concordance = external_mapped_feature_concordance(mapped, age)

        self.assertEqual(concordance.loc[0, "concordance"], "concordant")
        self.assertEqual(concordance.loc[0, "status"], "concordant_small_n_mapped")

    def test_external_module_contrasts_ignores_culture_samples(self):
        scores = pd.DataFrame(
            {
                "sample_class": ["biopsy", "biopsy", "culture"],
                "disease_group": ["healthy", "presbyosmia", "culture"],
                "module": ["m", "m", "m"],
                "mean_log1p_cpm": [1.0, 2.0, 10.0],
            }
        )

        contrasts = external_module_contrasts(scores)

        self.assertEqual(float(contrasts.loc[0, "presbyosmia_minus_healthy"]), 1.0)

    def test_external_validation_evidence_summary_claim_gates_small_n_results(self):
        config = {
            "datasets": {
                "oliva_2022": {
                    "title": "Olfactory aging",
                    "accession": "GSE184117",
                    "validation_use": "aging feature replication",
                    "expected_level": "single_cell",
                    "disease_context": ["aging", "presbyosmia"],
                    "source_url": "https://example.org/gse184117",
                },
                "bulk_context": {
                    "title": "Bulk context",
                    "validation_use": "marker sanity",
                    "expected_level": "bulk_rnaseq",
                    "disease_context": ["healthy_or_surgical_reference"],
                },
            }
        }
        dataset_summary = pd.DataFrame(
            {
                "dataset_id": ["oliva_2022", "bulk_context"],
                "readiness_class": ["ready_raw_adapter", "missing_files"],
                "ready_for_feature_validation": [False, False],
                "ready_for_raw_adapter": [True, False],
            }
        )
        metadata = pd.DataFrame(
            {
                "dataset_id": ["oliva_2022"] * 7,
                "sample_id": [f"GSM{i}" for i in range(7)],
                "donor_id": [f"D{i}" for i in range(7)],
                "disease_group": ["healthy", "healthy", "healthy", "presbyosmia", "presbyosmia", "presbyosmia", "culture"],
                "sample_class": ["biopsy", "biopsy", "biopsy", "biopsy", "biopsy", "biopsy", "culture"],
                "usable_for_external_validation": [True, True, True, True, True, True, False],
            }
        )
        module_contrasts = pd.DataFrame(
            {
                "module": ["regeneration"],
                "n_healthy": [3],
                "n_presbyosmia": [3],
                "status": ["descriptive_small_n"],
            }
        )
        marker_contrasts = pd.DataFrame(
            {
                "marker_panel": ["mature_osn"],
                "n_healthy": [3],
                "n_presbyosmia": [3],
                "status": ["marker_only_small_n"],
            }
        )

        evidence = external_validation_evidence_summary(
            config,
            dataset_summary,
            sample_metadata=metadata,
            module_contrasts=module_contrasts,
            marker_contrasts=marker_contrasts,
        )

        self.assertIn("raw_10x_sample_module_scores", evidence["evidence_type"].tolist())
        marker_row = evidence[evidence["evidence_type"].eq("raw_10x_marker_only_composition")].iloc[0]
        self.assertEqual(marker_row["supports_primary_claim"], "sanity_only")
        self.assertEqual(int(marker_row["n_samples"]), 6)
        bulk_row = evidence[evidence["dataset_id"].eq("bulk_context")].iloc[0]
        self.assertEqual(bulk_row["validation_strength"], "marker_context_only")

    def test_external_marker_age_concordance_maps_panels_to_gateway_features(self):
        marker_contrasts = pd.DataFrame(
            {
                "marker_panel": ["quiescent_hbc", "immune"],
                "n_healthy": [3, 3],
                "n_presbyosmia": [3, 3],
                "presbyosmia_minus_healthy": [0.2, -0.1],
                "p_value": [0.4, 0.8],
                "direction": ["higher_in_presbyosmia", "lower_in_presbyosmia"],
                "status": ["marker_only_small_n", "marker_only_small_n"],
            }
        )
        age_associations = pd.DataFrame(
            {
                "feature": ["prop__quiescent_hbc", "prop__inflammatory"],
                "direction": ["positive", "positive"],
                "beta_per_10_years": [0.1, 0.2],
                "p_value": [0.01, 0.02],
                "fdr": [0.05, 0.10],
            }
        )

        concordance = external_marker_age_concordance(marker_contrasts, age_associations)

        hbc = concordance[concordance["gateway_feature"].eq("prop__quiescent_hbc")].iloc[0]
        self.assertEqual(hbc["concordance"], "concordant")
        immune = concordance[concordance["gateway_feature"].eq("prop__inflammatory")].iloc[0]
        self.assertEqual(immune["concordance"], "discordant")
        self.assertIn("small_n_marker_only", hbc["status"])

def _write_gzip(path: Path, text: str) -> None:
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(text)


if __name__ == "__main__":
    unittest.main()
