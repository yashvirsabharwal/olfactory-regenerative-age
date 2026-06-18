#!/usr/bin/env python3
"""Validate a bounded scVI pilot output."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config
from ora.latent_recompute import validate_scvi_pilot
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/gateway.yaml")
    parser.add_argument("--h5ad", default=None)
    parser.add_argument("--embedding-key", default="X_scvi")
    parser.add_argument("--out", default=None)
    parser.add_argument("--max-validation-cells", type=int, default=50_000)
    parser.add_argument("--seed", type=int, default=13)
    args = parser.parse_args()

    config = load_config(args.config)
    outputs = config.get("outputs", {})
    h5ad = args.h5ad or outputs.get("latent_scvi_pilot_h5ad", "data/processed/gateway_scvi_pilot.h5ad")
    out = args.out or outputs.get("scvi_pilot_validation_tsv", "results/tables/scvi_pilot_validation.tsv")
    validation = validate_scvi_pilot(
        h5ad,
        embedding_key=args.embedding_key,
        max_validation_cells=args.max_validation_cells,
        seed=args.seed,
    )
    validation.to_csv(ensure_parent(out), sep="\t", index=False)
    print(f"Wrote scVI pilot validation: {out} ({validation.shape[0]} rows)")


if __name__ == "__main__":
    main()
