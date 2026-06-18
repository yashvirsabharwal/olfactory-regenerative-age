#!/usr/bin/env python3
"""Write a dry-run scVI latent recomputation feasibility plan."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config
from ora.latent_recompute import latent_recompute_feasibility, render_latent_recompute_workflow
from ora.reporting import load_schema
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/gateway.yaml")
    parser.add_argument("--schema", default=None)
    parser.add_argument("--h5ad", default=None)
    parser.add_argument("--output-h5ad", default=None)
    parser.add_argument("--n-top-genes", type=int, default=2000)
    parser.add_argument("--pilot-max-cells", type=int, default=25_000)
    parser.add_argument("--feasibility-out", default=None)
    parser.add_argument("--workflow-out", default=None)
    args = parser.parse_args()

    config = load_config(args.config)
    outputs = config.get("outputs", {})
    source = config.get("source", {})
    schema = load_schema(args.schema or outputs.get("schema_json", "results/reports/h5ad_schema.json"))
    h5ad_path = args.h5ad or source.get("h5ad_path", "data/raw/gateway.h5ad")
    output_h5ad = args.output_h5ad or outputs.get(
        "latent_scvi_pilot_h5ad",
        "data/processed/gateway_scvi_pilot.h5ad",
    )
    feasibility_out = args.feasibility_out or outputs.get(
        "latent_recompute_feasibility_tsv",
        "results/tables/latent_recompute_feasibility.tsv",
    )
    workflow_out = args.workflow_out or outputs.get(
        "latent_recompute_workflow_md",
        "docs/latent_recompute_workflow.md",
    )

    feasibility = latent_recompute_feasibility(
        schema,
        n_top_genes=args.n_top_genes,
        pilot_max_cells=args.pilot_max_cells,
    )
    workflow = render_latent_recompute_workflow(
        feasibility,
        h5ad_path=h5ad_path,
        output_h5ad=output_h5ad,
        n_top_genes=args.n_top_genes,
        pilot_max_cells=args.pilot_max_cells,
    )
    feasibility.to_csv(ensure_parent(feasibility_out), sep="\t", index=False)
    ensure_parent(workflow_out).write_text(workflow, encoding="utf-8")
    print(f"Wrote latent recompute feasibility: {feasibility_out} ({feasibility.shape[0]} rows)")
    print(f"Wrote latent recompute workflow: {workflow_out}")


if __name__ == "__main__":
    main()
