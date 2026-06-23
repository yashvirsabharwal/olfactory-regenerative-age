#!/usr/bin/env python3
"""Write a SOTA external validation candidate matrix from the registry."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config
from ora.external import external_candidate_matrix
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--external-config", default="configs/external_datasets.yaml")
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    config = load_config(args.external_config)
    outputs = config.get("outputs", {})
    out = args.out or outputs.get(
        "external_candidate_matrix_tsv",
        "results/tables/external_candidate_matrix.tsv",
    )
    matrix = external_candidate_matrix(config)
    matrix.to_csv(ensure_parent(out), sep="\t", index=False)
    print(f"Wrote external candidate matrix: {out} ({matrix.shape[0]} datasets)")


if __name__ == "__main__":
    main()
