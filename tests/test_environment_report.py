from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.environment import PackageSpec, lockfile_text, package_rows, r_environment_yml


class EnvironmentReportTests(unittest.TestCase):
    def test_package_rows_and_lockfile_mark_missing_versions(self):
        rows = package_rows(
            (PackageSpec("present-pkg", "present_module"), PackageSpec("missing-pkg")),
            group="toy",
            version_lookup=lambda package: "1.2.3" if package == "present-pkg" else None,
        )

        self.assertEqual(rows[0]["status"], "present")
        self.assertEqual(rows[1]["status"], "missing")
        text = lockfile_text(rows, title="Toy lock")
        self.assertIn("present-pkg==1.2.3", text)
        self.assertIn("# MISSING: missing-pkg", text)

    def test_r_environment_lists_required_bioconductor_packages(self):
        text = r_environment_yml()
        self.assertIn("bioconductor-edger", text)
        self.assertIn("bioconductor-limma", text)
        self.assertIn("bioconductor-milor", text)


if __name__ == "__main__":
    unittest.main()
