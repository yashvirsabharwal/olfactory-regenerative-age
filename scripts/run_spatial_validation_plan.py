#!/usr/bin/env python3
"""Build the spatial/histology validation design artifacts."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config
from ora.spatial_validation import (
    build_spatial_candidate_matrix,
    build_spatial_marker_panel,
    build_spatial_readout_plan,
    build_spatial_search_log,
    render_spatial_validation_plan,
    write_spatial_validation_outputs,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--external-config", default="configs/external_datasets.yaml")
    parser.add_argument("--marker-config", default="configs/spatial_validation_markers.yaml")
    parser.add_argument("--candidate-out", default="results/tables/spatial_validation_candidate_matrix.tsv")
    parser.add_argument("--marker-out", default="results/tables/spatial_validation_marker_panel.tsv")
    parser.add_argument("--readout-out", default="results/tables/spatial_validation_readout_plan.tsv")
    parser.add_argument("--search-log-out", default="results/tables/spatial_validation_search_log.tsv")
    parser.add_argument("--plan-out", default="docs/spatial_perturbation_validation_plan.md")
    args = parser.parse_args()

    external_config = load_config(args.external_config)
    marker_config = load_config(args.marker_config)
    candidate_matrix = build_spatial_candidate_matrix(external_config)
    marker_panel = build_spatial_marker_panel(marker_config)
    readout_plan = build_spatial_readout_plan(marker_config)
    search_log = build_spatial_search_log(external_config)
    plan_md = render_spatial_validation_plan(
        candidate_matrix=candidate_matrix,
        marker_panel=marker_panel,
        readout_plan=readout_plan,
        search_log=search_log,
    )
    write_spatial_validation_outputs(
        candidate_matrix=candidate_matrix,
        marker_panel=marker_panel,
        readout_plan=readout_plan,
        search_log=search_log,
        plan_md=plan_md,
        candidate_out=args.candidate_out,
        marker_out=args.marker_out,
        readout_out=args.readout_out,
        search_log_out=args.search_log_out,
        plan_out=args.plan_out,
    )
    usable = int(candidate_matrix["usable_for_primary_spatial_validation"].astype(bool).sum())
    print(
        "Wrote spatial validation plan: "
        f"{args.plan_out} ({candidate_matrix.shape[0]} candidates; {usable} primary-usable), "
        f"{args.marker_out} ({marker_panel.shape[0]} panels)"
    )


if __name__ == "__main__":
    main()
