from pathlib import Path
import sys
import unittest

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.metadata import parse_age_series


class AgeParsingTests(unittest.TestCase):
    def test_numeric_only_ages(self):
        parsed = parse_age_series(pd.Series([20, 35.5, None]))

        self.assertEqual(parsed.iloc[0], 20)
        self.assertEqual(parsed.iloc[1], 35.5)
        self.assertTrue(pd.isna(parsed.iloc[2]))

    def test_string_only_year_ages(self):
        parsed = parse_age_series(
            pd.Series(["65-year-old stage", "70 year", "85-years", "42 years old"])
        )

        self.assertEqual(parsed.tolist(), [65, 70, 85, 42])

    def test_mixed_numeric_string_and_missing_ages(self):
        parsed = parse_age_series(pd.Series([44, "65-year-old", "unknown", None, "72 yr old"]))

        self.assertEqual(parsed.iloc[0], 44)
        self.assertEqual(parsed.iloc[1], 65)
        self.assertTrue(pd.isna(parsed.iloc[2]))
        self.assertTrue(pd.isna(parsed.iloc[3]))
        self.assertEqual(parsed.iloc[4], 72)

    def test_malformed_values_do_not_parse(self):
        parsed = parse_age_series(pd.Series(["postnatal week 6", "adult", "age: sixty five", "2020 sample"]))

        self.assertTrue(parsed.isna().all())


if __name__ == "__main__":
    unittest.main()
