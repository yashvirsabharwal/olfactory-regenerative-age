#!/usr/bin/env python3
"""Inspect a Gateway H5AD in backed mode without loading the matrix."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config, project_path
from ora.io import inspect_h5ad
from ora.utils import ensure_parent, write_json


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/gateway.yaml")
    parser.add_argument("--h5ad", default=None)
    parser.add_argument("--schema-out", default=None)
    parser.add_argument("--obs-columns-out", default=None)
    parser.add_argument("--var-columns-out", default=None)
    args = parser.parse_args()

    config = load_config(args.config)
    h5ad_path = args.h5ad or config["source"]["h5ad_path"]
    schema, obs_table, var_table = inspect_h5ad(project_path(h5ad_path))

    outputs = config.get("outputs", {})
    schema_path = args.schema_out or outputs.get("schema_json", "results/reports/h5ad_schema.json")
    obs_columns_path = args.obs_columns_out or outputs.get("obs_columns_tsv", "data/metadata/gateway_obs_columns.tsv")
    var_columns_path = args.var_columns_out or outputs.get("var_columns_tsv", "data/metadata/gateway_var_columns.tsv")
    write_json(schema, schema_path)
    obs_table.to_csv(ensure_parent(obs_columns_path), sep="\t", index=False)
    var_table.to_csv(ensure_parent(var_columns_path), sep="\t", index=False)

    print(f"Inspected {h5ad_path}: {schema['n_obs']} cells x {schema['n_vars']} genes")
    print(f"Wrote {schema_path}")


if __name__ == "__main__":
    main()
