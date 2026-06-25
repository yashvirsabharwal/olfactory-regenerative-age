import unittest

import pandas as pd

from ora.organoid_perturbation import (
    gse309325_organoid_module_contrasts,
    parse_gse309325_sample_name,
)


class OrganoidPerturbationTests(unittest.TestCase):
    def test_parse_gse309325_sample_name(self):
        mock = parse_gse309325_sample_name("GSM9265687_hESC_1.csv.gz")
        infected = parse_gse309325_sample_name("GSM9265692_hESC_SARSCoV2_day7.csv.gz")

        self.assertEqual(mock["sample_id"], "GSM9265687")
        self.assertEqual(mock["condition"], "mock")
        self.assertEqual(mock["timepoint_day"], 0)
        self.assertEqual(mock["replicate"], 1)
        self.assertEqual(infected["condition"], "sars_cov_2")
        self.assertEqual(infected["timepoint_day"], 7)

    def test_organoid_module_contrasts_use_mock_baseline(self):
        scores = pd.DataFrame(
            {
                "sample_id": ["m1", "m2", "i1"],
                "sample_label": ["mock_1", "mock_2", "day1"],
                "condition": ["mock", "mock", "sars_cov_2"],
                "timepoint_day": [0, 0, 1],
                "module": ["ifn", "ifn", "ifn"],
                "mean_log1p_expression": [1.0, 3.0, 5.0],
            }
        )

        contrasts = gse309325_organoid_module_contrasts(scores)

        self.assertEqual(contrasts.shape[0], 1)
        self.assertEqual(float(contrasts.loc[0, "mock_mean"]), 2.0)
        self.assertEqual(float(contrasts.loc[0, "delta_vs_mock"]), 3.0)
        self.assertEqual(contrasts.loc[0, "direction_vs_mock"], "increased")


if __name__ == "__main__":
    unittest.main()
