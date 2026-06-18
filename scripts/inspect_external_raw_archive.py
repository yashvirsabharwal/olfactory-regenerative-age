#!/usr/bin/env python3
"""Inventory a raw external validation archive without extracting it."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config
from ora.external import inspect_external_archive
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gateway-config", default="configs/gateway.yaml")
    parser.add_argument("--archive", required=True)
    parser.add_argument("--dataset-id", default="external")
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    config = load_config(args.gateway_config)
    outputs = config.get("outputs", {})
    out_path = args.out or outputs.get("external_raw_inventory_tsv", "results/tables/external_raw_inventory.tsv")
    inventory = inspect_external_archive(args.archive, dataset_id=args.dataset_id)
    inventory.to_csv(ensure_parent(out_path), sep="\t", index=False)
    roles = ",".join(sorted(inventory["role"].dropna().astype(str).unique())) if not inventory.empty else "none"
    print(f"Wrote external raw inventory: {out_path} ({inventory.shape[0]} files; roles: {roles})")


if __name__ == "__main__":
    main()
