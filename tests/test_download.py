from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.download import infer_download_mode


class DownloadTests(unittest.TestCase):
    def test_download_mode_prefers_url(self):
        self.assertEqual(infer_download_mode("https://example.org/data.h5ad", "dataset"), "url")

    def test_download_mode_uses_dataset_id(self):
        self.assertEqual(infer_download_mode(None, "abc123"), "dataset_id")

    def test_download_mode_rejects_non_http_url(self):
        with self.assertRaises(ValueError):
            infer_download_mode("file:///tmp/data.h5ad", None)

    def test_download_mode_missing(self):
        self.assertEqual(infer_download_mode(None, None), "missing")


if __name__ == "__main__":
    unittest.main()

