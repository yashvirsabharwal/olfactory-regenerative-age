#!/usr/bin/env python3
"""Run frozen-score NDD label permutation within matched healthy context."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config
from ora.ndd import ndd_label_permutation
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gateway-config", default="configs/gateway.yaml")
    parser.add_argument("--projection", default=None)
    parser.add_argument("--n-permutations", type=int, default=5000)
    parser.add_argument("--random-seed", type=int, default=42)
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    config = load_config(args.gateway_config)
    outputs = config.get("outputs", {})
    projection_path = args.projection or outputs.get("ndd_ora_projection_tsv", "results/tables/ndd_ora_projection.tsv")
    out_path = args.out or outputs.get("ndd_label_permutation_tsv", "results/tables/ndd_label_permutation.tsv")
    projection = pd.read_csv(projection_path, sep="\t")
    result = ndd_label_permutation(
        projection,
        n_permutations=args.n_permutations,
        random_seed=args.random_seed,
    )
    result.to_csv(ensure_parent(out_path), sep="\t", index=False)
    print(f"Wrote NDD label permutation: {out_path} ({result.shape[0]} rows)")


if __name__ == "__main__":
    main()
