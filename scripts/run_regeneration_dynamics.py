#!/usr/bin/env python3
"""Run regeneration-dynamics feasibility audit and exploratory DPT pilot."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.regeneration_dynamics import (
    audit_h5ad_dynamics_inputs,
    build_dynamics_feasibility,
    render_regeneration_dynamics_report,
    run_scanpy_dpt_pilot,
    write_regeneration_dynamics_outputs,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--h5ad",
        action="append",
        default=[
            "data/raw/gateway.h5ad",
            "data/processed/gateway_scvi_lineage_basal_neural_100k.h5ad",
            "data/processed/foundation_benchmark_lineage_subset.h5ad",
            "data/processed/gateway_scvi_full_4m_reduced.h5ad",
        ],
        help="H5AD path to audit. Can be supplied multiple times.",
    )
    parser.add_argument(
        "--pilot-h5ad",
        default="data/processed/gateway_scvi_lineage_basal_neural_100k.h5ad",
    )
    parser.add_argument("--max-cells", type=int, default=40000)
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument("--audit-out", default="results/tables/regeneration_dynamics_input_audit.tsv")
    parser.add_argument(
        "--feasibility-out",
        default="results/tables/regeneration_dynamics_feasibility.tsv",
    )
    parser.add_argument("--report-out", default="results/reports/regeneration_dynamics_feasibility.md")
    parser.add_argument("--summary-out", default="results/tables/regeneration_dynamics_summary.tsv")
    parser.add_argument("--cells-out", default="results/tables/regeneration_dynamics_cell_scores.tsv")
    parser.add_argument("--figure-pdf", default="results/figures/extended_data_regeneration_dynamics.pdf")
    parser.add_argument("--figure-png", default="results/figures/extended_data_regeneration_dynamics.png")
    parser.add_argument("--skip-pilot", action="store_true")
    args = parser.parse_args()

    audit = audit_h5ad_dynamics_inputs(args.h5ad)
    feasibility = build_dynamics_feasibility(audit)
    dpt_summary = None
    dpt_cells = None
    if not args.skip_pilot:
        dpt_summary, dpt_cells = run_scanpy_dpt_pilot(
            h5ad_path=args.pilot_h5ad,
            max_cells=args.max_cells,
            seed=args.seed,
        )
    report = render_regeneration_dynamics_report(
        audit=audit,
        feasibility=feasibility,
        dpt_summary=dpt_summary,
    )
    write_regeneration_dynamics_outputs(
        audit=audit,
        feasibility=feasibility,
        report=report,
        audit_out=args.audit_out,
        feasibility_out=args.feasibility_out,
        report_out=args.report_out,
        dpt_summary=dpt_summary,
        dpt_cells=dpt_cells,
        summary_out=args.summary_out,
        cells_out=args.cells_out,
        figure_pdf=args.figure_pdf,
        figure_png=args.figure_png,
    )
    pilot_msg = ""
    if dpt_summary is not None:
        overall = dpt_summary[
            (dpt_summary["summary_type"] == "overall") & (dpt_summary["stratum"] == "all")
        ]
        if not overall.empty:
            pilot_msg = f"; DPT Spearman r={overall.iloc[0]['spearman_r']:.3f}"
    print(
        "Wrote regeneration dynamics feasibility: "
        f"{args.feasibility_out} ({feasibility.shape[0]} methods){pilot_msg}"
    )


if __name__ == "__main__":
    main()
