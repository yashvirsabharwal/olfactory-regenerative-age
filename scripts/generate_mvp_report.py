#!/usr/bin/env python3
"""Generate the composition-MVP ORA report and figures."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config
from ora.reporting import generate_mvp_report, load_schema


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/gateway.yaml")
    parser.add_argument("--manifest", default=None)
    parser.add_argument("--cohort-summary", default=None)
    parser.add_argument("--associations", default=None)
    parser.add_argument("--performance", default=None)
    parser.add_argument("--scores", default=None)
    parser.add_argument("--importance", default=None)
    parser.add_argument("--augmented-performance", default=None)
    parser.add_argument("--augmented-scores", default=None)
    parser.add_argument("--augmented-importance", default=None)
    parser.add_argument("--module-summary", default=None)
    parser.add_argument("--module-coverage", default=None)
    parser.add_argument("--donor-module-features", default=None)
    parser.add_argument("--schema", default=None)
    parser.add_argument("--out", default=None)
    parser.add_argument("--figure-dir", default=None)
    parser.add_argument("--top-n", type=int, default=12)
    args = parser.parse_args()

    config = load_config(args.config)
    outputs = config.get("outputs", {})
    paths = {
        "manifest": args.manifest or outputs.get("manifest_tsv", "data/processed/cohort_manifest.tsv"),
        "cohort_summary": args.cohort_summary or outputs.get("cohort_summary_tsv", "results/tables/cohort_summary.tsv"),
        "associations": args.associations or outputs.get("age_associations_tsv", "results/tables/age_cell_state_associations.tsv"),
        "performance": args.performance or outputs.get("model_performance_tsv", "results/tables/ora_model_performance.tsv"),
        "scores": args.scores or outputs.get("donor_ora_scores_tsv", "results/tables/donor_ora_scores.tsv"),
        "importance": args.importance or outputs.get("feature_importance_tsv", "results/tables/ora_feature_importance.tsv"),
        "augmented_performance": args.augmented_performance
        or outputs.get("augmented_model_performance_tsv", "results/tables/ora_augmented_model_performance.tsv"),
        "augmented_scores": args.augmented_scores
        or outputs.get("augmented_donor_ora_scores_tsv", "results/tables/augmented_donor_ora_scores.tsv"),
        "augmented_importance": args.augmented_importance
        or outputs.get("augmented_feature_importance_tsv", "results/tables/ora_augmented_feature_importance.tsv"),
        "module_summary": args.module_summary
        or outputs.get("module_score_summary_tsv", "results/tables/module_score_summary.tsv"),
        "module_coverage": args.module_coverage
        or outputs.get("module_gene_coverage_tsv", "results/tables/module_gene_coverage.tsv"),
        "donor_module_features": args.donor_module_features
        or outputs.get("donor_module_features_tsv", "data/processed/donor_module_features.tsv"),
        "schema": args.schema or outputs.get("schema_json", "results/reports/h5ad_schema.json"),
        "out": args.out or outputs.get("mvp_report_md", "results/reports/mvp_report.md"),
        "figure_dir": args.figure_dir or outputs.get("figure_dir", "results/figures"),
    }

    written = generate_mvp_report(
        manifest=pd.read_csv(paths["manifest"], sep="\t"),
        cohort_summary=pd.read_csv(paths["cohort_summary"], sep="\t"),
        associations=pd.read_csv(paths["associations"], sep="\t"),
        performance=pd.read_csv(paths["performance"], sep="\t"),
        scores=pd.read_csv(paths["scores"], sep="\t"),
        importance=pd.read_csv(paths["importance"], sep="\t"),
        augmented_performance=_read_optional_tsv(paths["augmented_performance"]),
        augmented_scores=_read_optional_tsv(paths["augmented_scores"]),
        augmented_importance=_read_optional_tsv(paths["augmented_importance"]),
        module_summary=_read_optional_tsv(paths["module_summary"]),
        module_coverage=_read_optional_tsv(paths["module_coverage"]),
        donor_module_features=_read_optional_tsv(paths["donor_module_features"]),
        schema=load_schema(paths["schema"]),
        source=config.get("source", {}),
        paper_defaults=config.get("paper_defaults", {}),
        out_md=paths["out"],
        figure_dir=paths["figure_dir"],
        top_n=args.top_n,
    )
    print(f"Wrote MVP report: {paths['out']}")
    print(f"Wrote {len(written) - 1} figures: {paths['figure_dir']}")


def _read_optional_tsv(path: str | None) -> pd.DataFrame | None:
    if not path:
        return None
    candidate = Path(path)
    if not candidate.exists():
        return None
    return pd.read_csv(candidate, sep="\t")


if __name__ == "__main__":
    main()
