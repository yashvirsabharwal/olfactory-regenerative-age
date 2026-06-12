#!/usr/bin/env python3
"""Summarize external validation readiness and published gene-list coverage."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config
from ora.external import external_dataset_summary, feature_matrix_contract_summary, published_gene_list_coverage
from ora.io import read_h5ad_backed
from ora.modules import DEFAULT_SYMBOL_COLUMNS
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--external-config", default="configs/external_datasets.yaml")
    parser.add_argument("--gateway-config", default="configs/gateway.yaml")
    parser.add_argument("--h5ad", default=None)
    parser.add_argument("--gene-table", default=None, help="Optional TSV/CSV with Gateway gene IDs or symbols if --h5ad is unavailable.")
    parser.add_argument("--summary-out", default=None)
    parser.add_argument("--coverage-out", default=None)
    parser.add_argument("--contract-out", default=None)
    args = parser.parse_args()

    external_config = load_config(args.external_config)
    gateway_config = load_config(args.gateway_config)
    outputs = external_config.get("outputs", {})
    summary_out = args.summary_out or outputs.get("validation_summary_tsv", "results/tables/external_validation_summary.tsv")
    coverage_out = args.coverage_out or outputs.get("gene_list_coverage_tsv", "results/tables/external_gene_list_coverage.tsv")
    contract_out = args.contract_out or outputs.get("feature_contract_tsv", "results/tables/external_feature_contract.tsv")

    dataset_summary = external_dataset_summary(external_config)
    dataset_summary.to_csv(ensure_parent(summary_out), sep="\t", index=False)
    feature_matrix_contract_summary(external_config).to_csv(ensure_parent(contract_out), sep="\t", index=False)

    var, var_names = _load_var(args.h5ad or gateway_config.get("source", {}).get("h5ad_path"), args.gene_table)
    symbol_columns = external_config.get("score", {}).get("var_symbol_columns", DEFAULT_SYMBOL_COLUMNS)
    coverage = published_gene_list_coverage(external_config, var, var_names, symbol_columns)
    coverage.to_csv(ensure_parent(coverage_out), sep="\t", index=False)

    ready = int(dataset_summary["ready_for_feature_validation"].sum()) if "ready_for_feature_validation" in dataset_summary else 0
    print(f"Wrote external validation summary: {summary_out} ({dataset_summary.shape[0]} datasets; {ready} feature-ready)")
    print(f"Wrote external feature contract: {contract_out}")
    print(f"Wrote published gene-list coverage: {coverage_out} ({coverage.shape[0]} gene lists)")


def _load_var(h5ad_path: str | None, gene_table_path: str | None) -> tuple[pd.DataFrame, pd.Index]:
    if h5ad_path:
        adata = read_h5ad_backed(h5ad_path)
        try:
            return adata.var.copy(), pd.Index(adata.var_names.astype(str))
        finally:
            close = getattr(adata, "file", None)
            if close is not None:
                close.close()
    if gene_table_path:
        sep = "," if str(gene_table_path).endswith(".csv") else "\t"
        table = pd.read_csv(gene_table_path, sep=sep)
        if "gene_id" in table:
            var = table.set_index("gene_id", drop=False)
            return var, pd.Index(var.index.astype(str))
        if "gene_symbol" in table:
            return table, pd.Index(table["gene_symbol"].astype(str))
        if "feature_name" in table:
            return table, pd.Index(table["feature_name"].astype(str))
        raise ValueError("--gene-table must include gene_id, gene_symbol, or feature_name.")
    raise ValueError("Provide --h5ad or --gene-table to resolve published gene-list coverage.")


if __name__ == "__main__":
    main()
