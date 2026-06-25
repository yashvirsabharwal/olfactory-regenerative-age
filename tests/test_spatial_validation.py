import unittest

from ora.spatial_validation import (
    build_spatial_candidate_matrix,
    build_spatial_marker_panel,
    build_spatial_search_log,
)


class SpatialValidationTests(unittest.TestCase):
    def test_candidate_matrix_adds_no_direct_spatial_sentinel_and_context_rows(self):
        config = {
            "public_data_exhaustion": {
                "candidates": [
                    {
                        "accession_or_dataset": "GSE235714",
                        "source_url": "https://example.org/gse235714",
                        "tissue": "Nasal tract / CRS tissue",
                        "assay": "NanoString GeoMx spatial",
                        "species": "human",
                        "age_availability": "not obvious",
                        "notes": "Specificity comparator.",
                    },
                    {
                        "accession_or_dataset": "GSE303809",
                        "tissue": "Fetal olfactory epithelium/head sections",
                        "assay": "MERFISH/spatial",
                        "species": "human fetal",
                        "age_availability": "developmental only; PCW 9 and 11",
                    },
                    {
                        "accession_or_dataset": "GSE184117",
                        "tissue": "Olfactory epithelium",
                        "assay": "10x scRNA-seq",
                        "species": "human",
                    },
                ]
            }
        }

        matrix = build_spatial_candidate_matrix(config)

        self.assertIn("direct_adult_human_olfactory_spatial_not_found", set(matrix["dataset_id"]))
        self.assertIn("gse235714", set(matrix["dataset_id"]))
        self.assertIn("gse303809", set(matrix["dataset_id"]))
        self.assertNotIn("gse184117", set(matrix["dataset_id"]))
        self.assertFalse(matrix["usable_for_primary_spatial_validation"].astype(bool).any())

    def test_marker_panel_flattens_markers(self):
        config = {
            "panels": [
                {
                    "panel_id": "hbc",
                    "theme": "basal",
                    "compartment": "basal_layer",
                    "priority": 1,
                    "expected_age_direction": "positive",
                    "markers": ["TP63", "KRT5"],
                    "readout": "density",
                    "rationale": "test",
                }
            ]
        }

        panel = build_spatial_marker_panel(config)

        self.assertEqual(panel.shape[0], 1)
        self.assertEqual(panel.loc[0, "markers"], "TP63,KRT5")
        self.assertEqual(int(panel.loc[0, "n_markers"]), 2)

    def test_search_log_keeps_spatial_queries_and_refresh_rows(self):
        config = {
            "public_data_exhaustion": {
                "search_log": [
                    {
                        "search_date": "2026-06-24",
                        "database_or_resource": "NCBI GEO",
                        "query_or_filter": "olfactory spatial transcriptomics",
                    },
                    {
                        "search_date": "2026-06-24",
                        "database_or_resource": "NCBI GEO",
                        "query_or_filter": "olfactory AND Homo sapiens",
                    },
                ]
            }
        }

        log = build_spatial_search_log(config)

        self.assertGreaterEqual(log.shape[0], 3)
        self.assertTrue(log["query_or_filter"].str.contains("spatial", case=False).any())
        self.assertTrue(log["search_date"].eq("2026-06-25").any())


if __name__ == "__main__":
    unittest.main()
