#!/usr/bin/env python3
"""Build a manuscript-oriented ORA feature interpretation table."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config
from ora.interpretation import build_feature_interpretation
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gateway-config", default="configs/gateway.yaml")
    parser.add_argument("--feature-stability", default=None)
    parser.add_argument("--associations", default=None)
    parser.add_argument("--out", default=None)
    parser.add_argument("--top-per-model", type=int, default=12)
    parser.add_argument("--top-n", type=int, default=30)
    args = parser.parse_args()

    config = load_config(args.gateway_config)
    outputs = config.get("outputs", {})
    stability_path = args.feature_stability or outputs.get(
        "ora_augmented_repeated_cv_feature_stability_tsv",
        "results/tables/ora_augmented_repeated_cv_feature_stability.tsv",
    )
    associations_path = args.associations or outputs.get(
        "age_associations_tsv",
        "results/tables/age_cell_state_associations.tsv",
    )
    out_path = args.out or outputs.get(
        "ora_feature_interpretation_tsv",
        "results/tables/ora_feature_interpretation.tsv",
    )

    interpretation = build_feature_interpretation(
        pd.read_csv(stability_path, sep="\t"),
        _read_optional_tsv(associations_path),
        top_per_model=args.top_per_model,
        top_n=args.top_n,
    )
    interpretation.to_csv(ensure_parent(out_path), sep="\t", index=False)
    print(f"Wrote ORA feature interpretation: {out_path} ({interpretation.shape[0]} rows)")


def _read_optional_tsv(path: str | None) -> pd.DataFrame | None:
    if not path:
        return None
    candidate = Path(path)
    if not candidate.exists():
        return None
    return pd.read_csv(candidate, sep="\t")


if __name__ == "__main__":
    main()
