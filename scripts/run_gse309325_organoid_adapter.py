#!/usr/bin/env python3
"""Score GSE309325 organoid perturbation modules."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config
from ora.organoid_perturbation import (
    render_gse309325_organoid_status,
    score_gse309325_organoid_modules,
    write_gse309325_organoid_outputs,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", default="data/external/GSE309325_RAW.tar")
    parser.add_argument("--gene-sets", default="configs/gene_sets.yaml")
    parser.add_argument("--regeneration-gene-sets", default="configs/regeneration_gene_sets.yaml")
    parser.add_argument("--chunksize", type=int, default=500)
    parser.add_argument("--qc-out", default="results/tables/gse309325_organoid_sample_qc.tsv")
    parser.add_argument("--scores-out", default="results/tables/gse309325_organoid_module_scores.tsv")
    parser.add_argument("--coverage-out", default="results/tables/gse309325_organoid_module_coverage.tsv")
    parser.add_argument("--contrasts-out", default="results/tables/gse309325_organoid_module_contrasts.tsv")
    parser.add_argument("--report-out", default="docs/gse309325_organoid_adapter_status.md")
    args = parser.parse_args()

    gene_set_config = _merged_gene_sets(args.gene_sets, args.regeneration_gene_sets)
    qc, scores, coverage, contrasts = score_gse309325_organoid_modules(
        args.archive,
        gene_set_config,
        chunksize=args.chunksize,
    )
    report = render_gse309325_organoid_status(
        qc=qc,
        scores=scores,
        coverage=coverage,
        contrasts=contrasts,
    )
    write_gse309325_organoid_outputs(
        qc=qc,
        scores=scores,
        coverage=coverage,
        contrasts=contrasts,
        report_md=report,
        qc_out=args.qc_out,
        scores_out=args.scores_out,
        coverage_out=args.coverage_out,
        contrasts_out=args.contrasts_out,
        report_out=args.report_out,
    )
    print(
        "Wrote GSE309325 organoid adapter outputs: "
        f"{args.contrasts_out} ({contrasts.shape[0]} contrasts), "
        f"{args.report_out}"
    )


def _merged_gene_sets(gene_sets_path: str, regeneration_gene_sets_path: str) -> dict:
    base = load_config(gene_sets_path)
    regen = load_config(regeneration_gene_sets_path)
    merged = {
        "score": {
            **base.get("score", {}),
            "log1p": True,
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
