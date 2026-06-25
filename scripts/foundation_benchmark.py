#!/usr/bin/env python3
"""Foundation-model benchmark planning and subset commands."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import project_path
from ora.foundation_benchmark import (
    DEFAULT_OUTPUTS,
    BenchmarkSubsetSpec,
    build_foundation_benchmark_subsets,
    load_configs,
)
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    _add_subsets(subparsers)
    _add_compare(subparsers)
    args = parser.parse_args()
    args.func(args)


def _add_subsets(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("subsets")
    parser.add_argument("--config", default="configs/gateway.yaml")
    parser.add_argument("--model-config", default="configs/models.yaml")
    parser.add_argument("--h5ad", default=None)
    parser.add_argument("--lineage-out", default=DEFAULT_OUTPUTS["lineage"])
    parser.add_argument("--epithelial-out", default=DEFAULT_OUTPUTS["epithelial"])
    parser.add_argument("--allcell-out", default=DEFAULT_OUTPUTS["allcell"])
    parser.add_argument("--manifest-out", default="results/tables/foundation_benchmark_subset_manifest.tsv")
    parser.add_argument("--donor-splits-out", default="results/tables/foundation_benchmark_donor_splits.tsv")
    parser.add_argument("--gene-manifest-out", default="results/tables/foundation_benchmark_gene_manifest.tsv")
    parser.add_argument("--gene-symbols-out", default="resources/foundation_benchmark/gateway_gene_symbols.txt")
    parser.add_argument("--gene-ids-out", default="resources/foundation_benchmark/gateway_gene_ids.txt")
    parser.add_argument("--lineage-cells", type=int, default=120_000)
    parser.add_argument("--epithelial-cells", type=int, default=180_000)
    parser.add_argument("--allcell-cells", type=int, default=250_000)
    parser.add_argument("--seed", type=int, default=20260625)
    parser.add_argument("--overwrite", action="store_true")
    parser.set_defaults(func=_subsets)


def _add_compare(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("compare")
    parser.add_argument("--geneformer-summary", default="results/tables/geneformer_age_model_summary.tsv")
    parser.add_argument("--runtime", default="results/tables/foundation_model_runtime.tsv")
    parser.add_argument("--out", default="results/tables/foundation_model_benchmark.tsv")
    parser.set_defaults(func=_compare)


def _subsets(args: argparse.Namespace) -> None:
    config, model_config = load_configs(args.config, args.model_config)
    h5ad = args.h5ad or config.get("source", {}).get("h5ad_path", "data/raw/gateway.h5ad")
    output_paths = {
        "lineage": args.lineage_out,
        "epithelial": args.epithelial_out,
        "allcell": args.allcell_out,
    }
    specs = (
        BenchmarkSubsetSpec(
            "lineage",
            args.lineage_cells,
            "Olfactory basal-to-neuronal lineage cells enriched for HBC, INP, iOSN, and mOSN states.",
        ),
        BenchmarkSubsetSpec(
            "epithelial",
            args.epithelial_cells,
            "Broad olfactory and respiratory epithelial subset, including sustentacular, secretory, glandular, and neuronal states.",
        ),
        BenchmarkSubsetSpec(
            "allcell",
            args.allcell_cells,
            "All-cell donor/fine-cell-type stratified subset for general-purpose foundation embeddings.",
        ),
    )
    manifest, splits, genes = build_foundation_benchmark_subsets(
        h5ad_path=project_path(h5ad),
        config=config,
        model_config=model_config,
        output_paths=output_paths,
        subset_specs=specs,
        manifest_out=args.manifest_out,
        donor_splits_out=args.donor_splits_out,
        gene_manifest_out=args.gene_manifest_out,
        gene_symbols_out=args.gene_symbols_out,
        gene_ids_out=args.gene_ids_out,
        seed=args.seed,
        overwrite=args.overwrite,
    )
    print(f"Wrote foundation benchmark subset manifest: {args.manifest_out} ({manifest.shape[0]} subsets)")
    print(f"Wrote donor split table: {args.donor_splits_out} ({splits.shape[0]} rows)")
    print(f"Wrote foundation benchmark gene manifest: {args.gene_manifest_out} ({genes.shape[0]} genes)")


def _compare(args: argparse.Namespace) -> None:
    rows: list[dict[str, object]] = []
    summary_path = Path(args.geneformer_summary)
    if summary_path.exists():
        summary = pd.read_csv(summary_path, sep="\t")
        if not summary.empty:
            best = summary.sort_values(["mae_mean", "rmse_mean"], ascending=True).iloc[0]
            rows.append(
                {
                    "model_family": "geneformer",
                    "checkpoint": "Geneformer-V1-10M",
                    "status": "ok",
                    "benchmark": "geneformer_v1_lineage_24k",
                    "model": best.get("model", ""),
                    "feature_set": "geneformer_v1_donor_embeddings",
                    "n": best.get("n", pd.NA),
                    "repeats": best.get("repeats", pd.NA),
                    "mae_mean": best.get("mae_mean", pd.NA),
                    "mae_ci_low": best.get("mae_ci_low", pd.NA),
                    "mae_ci_high": best.get("mae_ci_high", pd.NA),
                    "rmse_mean": best.get("rmse_mean", pd.NA),
                    "r2_mean": best.get("r2_mean", pd.NA),
                    "spearman_r_mean": best.get("spearman_r_mean", pd.NA),
                    "notes": "Best repeated donor-level age model using local Geneformer donor embeddings.",
                }
            )
    runtime_path = Path(args.runtime)
    if runtime_path.exists():
        runtime = pd.read_csv(runtime_path, sep="\t")
        for record in runtime.to_dict("records"):
            if record.get("status") == "ok":
                continue
            rows.append(
                {
                    "model_family": record.get("model_family", ""),
                    "checkpoint": record.get("checkpoint", ""),
                    "status": record.get("status", "not_run"),
                    "benchmark": "",
                    "model": "",
                    "feature_set": "",
                    "n": pd.NA,
                    "repeats": pd.NA,
                    "mae_mean": pd.NA,
                    "mae_ci_low": pd.NA,
                    "mae_ci_high": pd.NA,
                    "rmse_mean": pd.NA,
                    "r2_mean": pd.NA,
                    "spearman_r_mean": pd.NA,
                    "notes": record.get("notes", ""),
                }
            )
    table = pd.DataFrame(rows)
    table.to_csv(ensure_parent(args.out), sep="\t", index=False)
    print(f"Wrote foundation model benchmark table: {args.out} ({table.shape[0]} rows)")


if __name__ == "__main__":
    main()
