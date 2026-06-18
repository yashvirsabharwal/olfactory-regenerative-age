#!/usr/bin/env python3
"""Compare GSE184117 marker-only shifts with Gateway age-association directions."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config
from ora.external import external_marker_age_concordance
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gateway-config", default="configs/gateway.yaml")
    parser.add_argument("--marker-contrasts", default=None)
    parser.add_argument("--age-associations", default=None)
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    config = load_config(args.gateway_config)
    outputs = config.get("outputs", {})
    marker_path = args.marker_contrasts or outputs.get(
        "external_10x_marker_contrasts_tsv",
        "results/tables/external_10x_marker_contrasts.tsv",
    )
    age_path = args.age_associations or outputs.get(
        "age_associations_tsv",
        "results/tables/age_cell_state_associations.tsv",
    )
    out_path = args.out or outputs.get(
        "external_marker_age_concordance_tsv",
        "results/tables/external_marker_age_concordance.tsv",
    )
    concordance = external_marker_age_concordance(
        pd.read_csv(marker_path, sep="\t"),
        pd.read_csv(age_path, sep="\t"),
    )
    concordance.to_csv(ensure_parent(out_path), sep="\t", index=False)
    n_concordant = int(concordance["concordance"].eq("concordant").sum()) if not concordance.empty else 0
    n_discordant = int(concordance["concordance"].eq("discordant").sum()) if not concordance.empty else 0
    print(
        f"Wrote external marker-age concordance: {out_path} "
        f"({concordance.shape[0]} rows; {n_concordant} concordant, {n_discordant} discordant)"
    )


if __name__ == "__main__":
    main()
