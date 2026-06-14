#!/usr/bin/env python3
"""Summarize uncertainty and confounding context for NDD ORA projections."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config
from ora.ndd import summarize_ndd_projection_uncertainty
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gateway-config", default="configs/gateway.yaml")
    parser.add_argument("--projection", default=None)
    parser.add_argument("--n-bootstrap", type=int, default=5000)
    parser.add_argument("--random-seed", type=int, default=42)
    parser.add_argument("--uncertainty-out", default=None)
    parser.add_argument("--context-out", default=None)
    args = parser.parse_args()

    config = load_config(args.gateway_config)
    outputs = config.get("outputs", {})
    projection_path = args.projection or outputs.get("ndd_ora_projection_tsv", "results/tables/ndd_ora_projection.tsv")
    uncertainty_out = args.uncertainty_out or outputs.get(
        "ndd_ora_projection_uncertainty_tsv",
        "results/tables/ndd_ora_projection_uncertainty.tsv",
    )
    context_out = args.context_out or outputs.get(
        "ndd_ora_projection_context_tsv",
        "results/tables/ndd_ora_projection_context.tsv",
    )

    projection = pd.read_csv(projection_path, sep="\t")
    result = summarize_ndd_projection_uncertainty(
        projection,
        n_bootstrap=args.n_bootstrap,
        random_seed=args.random_seed,
    )
    result.uncertainty.to_csv(ensure_parent(uncertainty_out), sep="\t", index=False)
    result.context.to_csv(ensure_parent(context_out), sep="\t", index=False)
    print(f"Wrote NDD projection uncertainty: {uncertainty_out} ({result.uncertainty.shape[0]} rows)")
    print(f"Wrote NDD projection context: {context_out} ({result.context.shape[0]} rows)")


if __name__ == "__main__":
    main()
