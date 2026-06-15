#!/usr/bin/env python3
"""Summarize AD/PD ORA projection diagnostics by covariate and yield strata."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config
from ora.ndd import ndd_projection_diagnostics
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gateway-config", default="configs/gateway.yaml")
    parser.add_argument("--projection", default=None)
    parser.add_argument("--out", default=None)
    parser.add_argument("--min-donors-ok", type=int, default=2)
    args = parser.parse_args()

    config = load_config(args.gateway_config)
    outputs = config.get("outputs", {})
    projection_path = args.projection or outputs.get("ndd_ora_projection_tsv", "results/tables/ndd_ora_projection.tsv")
    out_path = args.out or outputs.get(
        "ndd_ora_projection_diagnostics_tsv",
        "results/tables/ndd_ora_projection_diagnostics.tsv",
    )
    projection = pd.read_csv(projection_path, sep="\t")
    diagnostics = ndd_projection_diagnostics(projection, min_donors_ok=args.min_donors_ok)
    diagnostics.to_csv(ensure_parent(out_path), sep="\t", index=False)
    print(f"Wrote NDD projection diagnostics: {out_path} ({diagnostics.shape[0]} rows)")


if __name__ == "__main__":
    main()
