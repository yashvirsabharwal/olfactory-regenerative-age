#!/usr/bin/env python3
"""Run donor-level cross-tissue CELLxGENE age-effect estimates."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config
from ora.cross_tissue_age import (
    build_cellxgene_asset_inventory,
    build_cellxgene_donor_feature_matrix,
    estimate_cross_tissue_age_effects,
    parse_dataset_specs,
    render_cross_tissue_age_effect_report,
    summarize_ora_cross_tissue_age_effects,
    write_cross_tissue_age_outputs,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/cross_tissue_age_effects.yaml")
    parser.add_argument("--gene-sets", default="configs/gene_sets.yaml")
    parser.add_argument("--regeneration-gene-sets", default="configs/regeneration_gene_sets.yaml")
    parser.add_argument(
        "--classification",
        default="results/tables/ora_cross_tissue_feature_classification.tsv",
    )
    parser.add_argument(
        "--inventory-out",
        default="results/tables/cross_tissue_cellxgene_asset_inventory.tsv",
    )
    parser.add_argument(
        "--donor-features-out",
        default="data/processed/cross_tissue_cellxgene_donor_features.tsv",
    )
    parser.add_argument(
        "--module-coverage-out",
        default="results/tables/cross_tissue_cellxgene_module_coverage.tsv",
    )
    parser.add_argument(
        "--effects-out",
        default="results/tables/cross_tissue_age_effects.tsv",
    )
    parser.add_argument(
        "--ora-summary-out",
        default="results/tables/ora_cross_tissue_age_effect_summary.tsv",
    )
    parser.add_argument("--report-out", default="docs/cross_tissue_age_effects.md")
    args = parser.parse_args()

    config = load_config(args.config)
    settings = config.get("settings", {})
    gene_set_config = _merged_gene_sets(args.gene_sets, args.regeneration_gene_sets)
    dataset_specs = parse_dataset_specs(config)
    if not dataset_specs:
        raise SystemExit("No cross-tissue datasets were configured.")

    inventory = build_cellxgene_asset_inventory(dataset_specs)
    donor_features, module_coverage = build_cellxgene_donor_feature_matrix(
        dataset_specs,
        gene_set_config,
        chunk_size=int(settings.get("chunk_size", 20_000)),
    )
    effects = estimate_cross_tissue_age_effects(
        donor_features,
        min_donors=int(settings.get("min_donors", 4)),
        min_cells_per_donor=int(settings.get("min_cells_per_donor", 20)),
        adult_min_age=float(settings.get("adult_min_age_years", 18)),
    )
    classification = pd.read_csv(args.classification, sep="\t") if Path(args.classification).exists() else pd.DataFrame()
    ora_summary = summarize_ora_cross_tissue_age_effects(classification, effects)
    report = render_cross_tissue_age_effect_report(
        inventory=inventory,
        effects=effects,
        ora_summary=ora_summary,
    )
    write_cross_tissue_age_outputs(
        inventory=inventory,
        donor_features=donor_features,
        module_coverage=module_coverage,
        effects=effects,
        ora_summary=ora_summary,
        report_md=report,
        inventory_out=args.inventory_out,
        donor_features_out=args.donor_features_out,
        module_coverage_out=args.module_coverage_out,
        effects_out=args.effects_out,
        ora_summary_out=args.ora_summary_out,
        report_out=args.report_out,
    )
    adult_ok = int(
        effects["analysis_scope"].eq("adult_only").mul(effects["status"].eq("ok")).sum()
        if not effects.empty
        else 0
    )
    print(
        "Wrote cross-tissue age effects: "
        f"{args.effects_out} ({effects.shape[0]} rows; {adult_ok} adult estimates), "
        f"{args.ora_summary_out} ({ora_summary.shape[0]} ORA features)"
    )


def _merged_gene_sets(gene_sets_path: str, regeneration_gene_sets_path: str) -> dict:
    base = load_config(gene_sets_path)
    regen = load_config(regeneration_gene_sets_path)
    merged = {
        "score": {
            **base.get("score", {}),
            "log1p": False,
        },
        "gene_sets": {},
    }
    merged["gene_sets"].update(base.get("gene_sets", {}))
    for name, spec in regen.get("gene_sets", {}).items():
        key = str(name)
        if key in merged["gene_sets"]:
            key = f"regeneration_{key}"
        merged["gene_sets"][key] = spec
    return merged


if __name__ == "__main__":
    main()
