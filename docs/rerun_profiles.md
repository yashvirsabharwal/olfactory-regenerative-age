# Rerun Profiles

Updated: 2026-06-24

This document defines the staged rerun profiles used for review, collaboration, and submission freeze. The Makefile is the canonical workflow entry point.

## Profile Summary

| Profile | Command | Use case | Expected runtime | Compute profile |
| --- | --- | --- | --- | --- |
| `local-light` | `PYTHON=.venv/bin/python make local-light` | Fast reviewer/collaborator check of tests, lint, external summaries, manuscript tables/figures, manuscript package checks, and provenance. | Minutes on the M5 laptop if result artifacts already exist. | CPU, modest memory. |
| `local-medium` | `PYTHON=.venv/bin/python make local-medium` | Regenerate core donor-level ORA summaries, compositional/negative-control analyses, model diagnostics, external evidence ledger, NDD guardrails, then run `local-light`. | Tens of minutes to hours depending on model runs. | CPU, 32 GB RAM laptop acceptable; optional native OpenMP helps XGBoost/LightGBM. |
| `remote-heavy` | `PYTHON=.venv/bin/python make remote-heavy` | Recompute full-scale latent and neighborhood claims, including all-cell scVI, Milo-style full/matched lineage analyses, edgeR parity, and MiloR subset sensitivity. | Hours to days. | Remote CPU/GPU strongly recommended; large memory and stable storage required. |
| `submission-freeze` | `PYTHON=.venv/bin/python make submission-freeze` | Final non-TeX freeze: environment report, release manifest, archive review package, external status reports, publication tables, figures, manuscript package check, and provenance. | Minutes if artifacts exist. | CPU, modest memory. |

## `local-light`

Purpose: verify that the repository and already-generated manuscript-facing artifacts are internally consistent.

Expected inputs:

- Source code, configs, and tests.
- Existing result tables under `results/tables`.
- Existing manuscript figures or inputs needed by `scripts/build_manuscript_figures.py`.
- Existing large generated objects referenced by provenance; deferred heavy artifacts may remain deferred.

Primary outputs:

- Unit-test and Ruff pass/fail status.
- `results/tables/external_candidate_matrix.tsv`
- `results/tables/public_data_exhaustion.tsv`
- `results/reports/public_data_exhaustion.md`
- `results/reports/gse184117_reanalysis_status.md`
- `results/reports/publication_tables.md`
- `results/reports/manuscript_package_check.md`
- `results/reports/output_provenance.tsv`

Stop conditions:

- Any test or lint failure.
- Manuscript package check reports missing non-deferred figure/table/citation assets.
- Output provenance reports missing non-deferred outputs.

## `local-medium`

Purpose: refresh the donor-level statistical core without recomputing the all-cell latent atlas.

Expected inputs:

- `data/raw/gateway.h5ad`
- `data/external/GSE184117_RAW.tar`
- `data/external/GSE184117_series_matrix.txt.gz`
- Existing model configs in `configs/models.yaml`
- Existing Gateway and external configs.

Primary outputs:

- `data/processed/cohort_manifest.tsv`
- `data/processed/donor_cell_state_features.tsv`
- `data/processed/ora_feature_matrix.tsv`
- `data/processed/ora_augmented_feature_matrix.tsv`
- ORA model diagnostics, repeated CV, permutation-null, tuning, stacking, and feature-interpretation tables.
- Compositional and negative-control tables/figures.
- NDD projection guardrail tables.
- Current external evidence ledger.
- All `local-light` outputs.

Stop conditions:

- Backend fallback occurs without an explicit `--allow-fallback` path in the relevant target.
- Donor/sample counts deviate from expected cohort summary without a documented source-data change.
- Negative-control or technical-baseline outputs change the manuscript claim boundary.

## `remote-heavy`

Purpose: refresh the expensive latent/neighborhood layer that supports the secondary all-cell spatial-neighborhood claim.

Expected inputs:

- `data/raw/gateway.h5ad`
- `resources/scvi/full_4m_genes.txt`
- Existing cohort manifest.
- R/Bioconductor/MiloR environment configured through the documented environment files.
- Enough disk space for reduced and mapped H5AD files plus model directories.

Primary outputs:

- `data/processed/gateway_hvg3003_4m.h5ad`
- `data/processed/gateway_scvi_full_4m_reduced.h5ad`
- `results/models/gateway_scvi_full_4m_reduced/`
- `results/tables/scvi_full_4m_reduced_validation.tsv`
- `results/tables/scvi_embedding_claim_gates.tsv`
- Full and matched lineage Milo-style neighborhood, program, age-bin, edgeR parity, and MiloR subset sensitivity tables.

Resource notes:

- Use remote CPU/GPU for all-cell scVI and full-neighborhood work.
- The laptop can prepare commands and validate small outputs, but stable remote storage is preferred for the full 4M artifacts.
- Before running, confirm the archive/reviewer-access destination so large outputs can be staged once regenerated.

Stop conditions:

- Full scVI training fails or produces missing/invalid latent embeddings.
- Milo-style and edgeR parity outputs disagree in a way that changes the latent-neighborhood claim.
- Official MiloR subset sensitivity cannot run and is not documented as blocked with environment details.

## `submission-freeze`

Purpose: produce the final non-TeX review package after analysis artifacts are frozen.

Expected inputs:

- Current generated tables, figures, model outputs, and H5AD/model artifacts.
- Final external evidence and public-data exhaustion reports.
- Current command manifest and release artifact registry.
- Confirmed archive/reviewer-access policy, if available.

Primary outputs:

- `results/reports/environment_report.md`
- `results/reports/release_artifact_manifest.md`
- `results/reports/archive_review_package.md`
- `results/reports/gse184117_reanalysis_status.md`
- `results/reports/public_data_exhaustion.md`
- `results/reports/publication_tables.md`
- `results/reports/manuscript_package_check.md`
- `results/reports/output_provenance.tsv`

Stop conditions:

- Archive review package still reports blockers that are not acceptable for the target journal.
- Manuscript package check reports unresolved assets or citations.
- Provenance reports missing non-deferred outputs.
- Local TeX remains unavailable when a compiled PDF is required.

## TeX PDF Build

`submission-freeze` deliberately stops before `make manuscript` because the laptop currently lacks a TeX engine. Build the PDF in a TeX-enabled environment with:

```bash
PYTHON=.venv/bin/python make manuscript
```

The PDF build is required before final submission and remains a Gate A blocker until `manuscript/main.pdf` is compiled and visually inspected.
