#!/usr/bin/env python3
"""Download the Gateway H5AD from CELLxGENE or a direct source URL."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config, project_path
from ora.download import download_cellxgene_dataset, download_direct_url, infer_download_mode


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/gateway.yaml")
    parser.add_argument("--dataset-id", default=None)
    parser.add_argument("--url", default=None, help="Direct H5AD download URL.")
    parser.add_argument("--out", default=None)
    parser.add_argument("--dry-run", action="store_true", help="Print resolved download source without downloading.")
    parser.add_argument("--force", action="store_true", help="Overwrite an existing output file.")
    args = parser.parse_args()

    config = load_config(args.config)
    source = config.get("source", {})
    dataset_id = args.dataset_id or source.get("dataset_id")
    output_path = project_path(args.out or source.get("h5ad_path", "data/raw/gateway.h5ad"))
    direct_url = args.url or source.get("direct_h5ad_url")
    mode = infer_download_mode(direct_url, dataset_id)

    if mode == "missing":
        collection_id = source.get("collection_id")
        raise SystemExit(
            "No dataset ID or direct H5AD URL was provided. "
            f"Configured collection_id is {collection_id!r}; resolve it to a concrete dataset ID "
            "or pass --url from CELLxGENE, then rerun this command."
        )
    expected = source.get("h5ad_filesize_bytes")
    if args.dry_run:
        print(f"mode: {mode}")
        print(f"dataset_id: {dataset_id}")
        print(f"url: {direct_url}")
        print(f"output: {output_path}")
        if expected:
            print(f"expected_bytes: {expected}")
        return
    if output_path.exists() and not args.force:
        raise SystemExit(f"Output already exists: {output_path}. Pass --force to overwrite.")
    if mode == "url":
        path = download_direct_url(str(direct_url), output_path)
    else:
        path = download_cellxgene_dataset(str(dataset_id), output_path)
    print(f"Wrote H5AD: {path}")


if __name__ == "__main__":
    main()
