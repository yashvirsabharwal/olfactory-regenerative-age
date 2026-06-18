#!/usr/bin/env python3
"""Write command and output provenance tables for generated ORA artifacts."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config
from ora.provenance import command_manifest_table, output_provenance_table
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gateway-config", default="configs/gateway.yaml")
    parser.add_argument("--command-manifest", default="configs/command_manifest.yaml")
    parser.add_argument("--checksum-max-bytes", type=int, default=100 * 1024 * 1024)
    parser.add_argument("--command-out", default=None)
    parser.add_argument("--provenance-out", default=None)
    args = parser.parse_args()

    gateway_config = load_config(args.gateway_config)
    command_config = load_config(args.command_manifest)
    outputs = gateway_config.get("outputs", {})
    command_out = args.command_out or outputs.get("command_manifest_tsv", "results/reports/command_manifest.tsv")
    provenance_out = args.provenance_out or outputs.get("output_provenance_tsv", "results/reports/output_provenance.tsv")
    commands = command_manifest_table(command_config)
    provenance = output_provenance_table(command_config, checksum_max_bytes=args.checksum_max_bytes)
    commands.to_csv(ensure_parent(command_out), sep="\t", index=False)
    provenance.to_csv(ensure_parent(provenance_out), sep="\t", index=False)
    missing = int((~provenance["exists"].astype(bool)).sum()) if "exists" in provenance else 0
    print(f"Wrote command manifest: {command_out} ({commands.shape[0]} stages)")
    print(f"Wrote output provenance: {provenance_out} ({provenance.shape[0]} outputs; {missing} missing)")


if __name__ == "__main__":
    main()
