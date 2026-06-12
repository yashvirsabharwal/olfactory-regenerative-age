#!/usr/bin/env python3
"""Run the MVP pipeline on a temporary toy Gateway-shaped H5AD."""

from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    python = sys.executable
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        h5ad = tmp / "toy_gateway.h5ad"
        schema = tmp / "h5ad_schema.json"
        obs_columns = tmp / "obs_columns.tsv"
        var_columns = tmp / "var_columns.tsv"
        manifest = tmp / "cohort_manifest.tsv"
        cohort_summary = tmp / "cohort_summary.tsv"
        counts = tmp / "counts.tsv"
        cell_features = tmp / "cell_features.tsv"
        matrix = tmp / "ora_feature_matrix.tsv"
        augmented_matrix = tmp / "ora_augmented_feature_matrix.tsv"
        associations = tmp / "age_associations.tsv"
        performance = tmp / "ora_performance.tsv"
        scores = tmp / "ora_scores.tsv"
        importance = tmp / "importance.tsv"
        augmented_performance = tmp / "ora_augmented_performance.tsv"
        augmented_scores = tmp / "augmented_ora_scores.tsv"
        augmented_importance = tmp / "augmented_importance.tsv"
        ndd_projection = tmp / "ndd_ora_projection.tsv"
        ndd_projection_summary = tmp / "ndd_ora_projection_summary.tsv"
        module_summary = tmp / "module_score_summary.tsv"
        module_coverage = tmp / "module_gene_coverage.tsv"
        donor_module_features = tmp / "donor_module_features.tsv"
        pseudobulk_counts = tmp / "pseudobulk_counts.tsv.gz"
        pseudobulk_metadata = tmp / "pseudobulk_metadata.tsv"
        pseudobulk_coverage = tmp / "pseudobulk_gene_coverage.tsv"
        pseudobulk_de = tmp / "pseudobulk_de.tsv"
        pseudobulk_covariate_de = tmp / "pseudobulk_covariate_de.tsv"
        report = tmp / "mvp_report.md"
        figure_dir = tmp / "figures"

        _run([python, "scripts/create_toy_gateway_h5ad.py", "--out", str(h5ad)], root)
        _run(
            [
                python,
                "scripts/inspect_h5ad.py",
                "--h5ad",
                str(h5ad),
                "--schema-out",
                str(schema),
                "--obs-columns-out",
                str(obs_columns),
                "--var-columns-out",
                str(var_columns),
            ],
            root,
        )
        _run(
            [
                python,
                "scripts/build_sample_manifest.py",
                "--h5ad",
                str(h5ad),
                "--out",
                str(manifest),
                "--summary-out",
                str(cohort_summary),
            ],
            root,
        )
        _run(
            [
                python,
                "scripts/aggregate_cell_counts.py",
                "--h5ad",
                str(h5ad),
                "--counts-out",
                str(counts),
                "--features-out",
                str(cell_features),
            ],
            root,
        )
        _run([python, "scripts/build_feature_matrix.py", "--features", str(cell_features), "--out", str(matrix)], root)
        _run(
            [
                python,
                "scripts/run_age_associations.py",
                "--features",
                str(cell_features),
                "--manifest",
                str(manifest),
                "--out",
                str(associations),
            ],
            root,
        )
        _run(
            [
                python,
                "scripts/run_age_models.py",
                "--features",
                str(matrix),
                "--manifest",
                str(manifest),
                "--performance-out",
                str(performance),
                "--scores-out",
                str(scores),
                "--importance-out",
                str(importance),
            ],
            root,
        )
        _run(
            [
                python,
                "scripts/score_gene_sets.py",
                "--h5ad",
                str(h5ad),
                "--summary-out",
                str(module_summary),
                "--coverage-out",
                str(module_coverage),
                "--donor-features-out",
                str(donor_module_features),
                "--chunk-size",
                "10",
            ],
            root,
        )
        _run(
            [
                python,
                "scripts/aggregate_pseudobulk.py",
                "--h5ad",
                str(h5ad),
                "--counts-out",
                str(pseudobulk_counts),
                "--metadata-out",
                str(pseudobulk_metadata),
                "--coverage-out",
                str(pseudobulk_coverage),
                "--de-out",
                str(pseudobulk_de),
                "--chunk-size",
                "10",
                "--min-donors",
                "2",
            ],
            root,
        )
        _run(
            [
                python,
                "scripts/run_pseudobulk_covariate_de.py",
                "--counts",
                str(pseudobulk_counts),
                "--metadata",
                str(pseudobulk_metadata),
                "--coverage",
                str(pseudobulk_coverage),
                "--manifest",
                str(manifest),
                "--out",
                str(pseudobulk_covariate_de),
                "--min-donors",
                "2",
            ],
            root,
        )
        _run(
            [
                python,
                "scripts/build_feature_matrix.py",
                "--features",
                str(cell_features),
                "--module-features",
                str(donor_module_features),
                "--out",
                str(augmented_matrix),
            ],
            root,
        )
        _run(
            [
                python,
                "scripts/run_age_models.py",
                "--features",
                str(augmented_matrix),
                "--manifest",
                str(manifest),
                "--performance-out",
                str(augmented_performance),
                "--scores-out",
                str(augmented_scores),
                "--importance-out",
                str(augmented_importance),
            ],
            root,
        )
        _run(
            [
                python,
                "scripts/project_ndd_ora.py",
                "--features",
                str(augmented_matrix),
                "--manifest",
                str(manifest),
                "--scores-out",
                str(ndd_projection),
                "--summary-out",
                str(ndd_projection_summary),
            ],
            root,
        )
        _run(
            [
                python,
                "scripts/generate_mvp_report.py",
                "--manifest",
                str(manifest),
                "--cohort-summary",
                str(cohort_summary),
                "--associations",
                str(associations),
                "--performance",
                str(performance),
                "--scores",
                str(scores),
                "--importance",
                str(importance),
                "--augmented-performance",
                str(augmented_performance),
                "--augmented-scores",
                str(augmented_scores),
                "--augmented-importance",
                str(augmented_importance),
                "--ndd-projection",
                str(ndd_projection),
                "--ndd-projection-summary",
                str(ndd_projection_summary),
                "--module-summary",
                str(module_summary),
                "--module-coverage",
                str(module_coverage),
                "--donor-module-features",
                str(donor_module_features),
                "--pseudobulk-de",
                str(pseudobulk_de),
                "--pseudobulk-coverage",
                str(pseudobulk_coverage),
                "--pseudobulk-metadata",
                str(pseudobulk_metadata),
                "--pseudobulk-covariate-de",
                str(pseudobulk_covariate_de),
                "--schema",
                str(schema),
                "--out",
                str(report),
                "--figure-dir",
                str(figure_dir),
            ],
            root,
        )
        print(f"Toy smoke workflow completed in {tmp}")


def _run(command: list[str], cwd: Path) -> None:
    print("+ " + " ".join(command))
    subprocess.run(command, cwd=cwd, check=True)


if __name__ == "__main__":
    main()
