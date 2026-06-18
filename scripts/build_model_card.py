#!/usr/bin/env python3
"""Build a model-card table from existing ORA benchmark summaries."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config
from ora.manuscript import build_model_card
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gateway-config", default="configs/gateway.yaml")
    parser.add_argument("--feature-set-comparison", default=None)
    parser.add_argument("--calibration", default=None)
    parser.add_argument("--permutation", default=None)
    parser.add_argument("--nested-tuning", default=None)
    parser.add_argument("--stacking", default=None)
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    config = load_config(args.gateway_config)
    outputs = config.get("outputs", {})
    paths = {
        "feature_set_comparison": args.feature_set_comparison
        or outputs.get("ora_feature_set_model_comparison_tsv", "results/tables/ora_feature_set_model_comparison.tsv"),
        "calibration": args.calibration or outputs.get("ora_calibration_tsv", "results/tables/ora_calibration.tsv"),
        "permutation": args.permutation or outputs.get("ora_permutation_empirical_tsv", "results/tables/ora_permutation_empirical.tsv"),
        "nested_tuning": args.nested_tuning or outputs.get("ora_nested_tuning_summary_tsv", "results/tables/ora_nested_tuning_summary.tsv"),
        "stacking": args.stacking or outputs.get("ora_stacking_summary_tsv", "results/tables/ora_stacking_summary.tsv"),
    }
    card = build_model_card(
        feature_set_comparison=_read_optional(paths["feature_set_comparison"]),
        calibration=_read_optional(paths["calibration"]),
        permutation=_read_optional(paths["permutation"]),
        nested_tuning=_read_optional(paths["nested_tuning"]),
        stacking=_read_optional(paths["stacking"]),
    )
    out_path = args.out or outputs.get("ora_model_card_tsv", "results/tables/ora_model_card.tsv")
    card.to_csv(ensure_parent(out_path), sep="\t", index=False)
    print(f"Wrote ORA model card: {out_path} ({card.shape[0]} rows)")


def _read_optional(path: str | None) -> pd.DataFrame | None:
    if not path:
        return None
    candidate = Path(path)
    if not candidate.exists():
        return None
    return pd.read_csv(candidate, sep="\t")


if __name__ == "__main__":
    main()
