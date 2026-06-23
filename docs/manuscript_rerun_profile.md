# Manuscript Rerun Profile

Updated: 2026-06-23

Purpose: provide a staged, reviewer-readable rerun order for manuscript-facing artifacts. The Makefile is the canonical workflow; this profile describes the recommended order and separates local-ready work from heavyweight deferred runs.

## Profile A: Local Manuscript Refresh

Use this profile when raw and processed inputs already exist locally and the goal is to refresh manuscript-facing tables, figures, and provenance without retraining the full 4M latent substrate.

| Stage | Command | Expected role |
| --- | --- | --- |
| Environment smoke check | `.venv/bin/python -m ruff check .` and `PYTHON=.venv/bin/python make test` | Confirms code and tests before regenerating outputs. |
| External registry/evidence refresh | `PYTHON=.venv/bin/python make external-validation` then `PYTHON=.venv/bin/python make external-evidence` | Rebuilds external readiness and claim-strength ledgers from `configs/external_datasets.yaml`. |
| ORA summary tables | `PYTHON=.venv/bin/python make model-card` and `PYTHON=.venv/bin/python make feature-interpretation` | Rebuilds compact model-card and interpretation tables. |
| Guardrail summaries | `PYTHON=.venv/bin/python make project-ndd-uncertainty`, `PYTHON=.venv/bin/python make project-ndd-diagnostics`, `PYTHON=.venv/bin/python make project-ndd-label-permutation` | Refreshes exploratory AD/PD projection guardrails. |
| Publication tables | `PYTHON=.venv/bin/python make publication-tables` | Rebuilds `results/tables/manuscript_table_*.tsv` and `results/reports/publication_tables.md`. |
| Manuscript figures | `PYTHON=.venv/bin/python make manuscript-figures` | Rebuilds main and extended-data figure files from current result tables. |
| Provenance | `PYTHON=.venv/bin/python make output-provenance` | Rebuilds command manifest and output provenance. |
| Report | `PYTHON=.venv/bin/python make` | Optional local Markdown refresh. |

## Profile B: Heavyweight Latent/Neighborhood Refresh

Use this profile on a remote or high-memory workstation when the 4M latent substrate must be restored or regenerated before final submission.

| Stage | Command | Required artifact |
| --- | --- | --- |
| Reduced 4M substrate | `PYTHON=.venv/bin/python make scvi-reduced-4m` | `data/processed/gateway_hvg3003_4m.h5ad` |
| Primary full 4M reduced scVI | `PYTHON=.venv/bin/python make scvi-full-4m-reduced` | `data/processed/gateway_scvi_full_4m_reduced.h5ad`, `results/models/gateway_scvi_full_4m_reduced` |
| Full 4M validation | `PYTHON=.venv/bin/python make scvi-full-validation` | `results/tables/scvi_full_4m_reduced_validation.tsv` |
| Full 4M lineage Milo-style | `PYTHON=.venv/bin/python make milo-full-4m-lineage` and `PYTHON=.venv/bin/python make milo-full-4m-lineage-matched` | Full and matched lineage neighborhood DA summaries. |
| Program and age-bin summaries | `PYTHON=.venv/bin/python make milo-full-4m-lineage-programs`, `PYTHON=.venv/bin/python make milo-full-4m-lineage-matched-programs`, `PYTHON=.venv/bin/python make milo-full-4m-lineage-age-bins`, `PYTHON=.venv/bin/python make milo-full-4m-lineage-matched-age-bins` | Program enrichment and age-bin robustness tables. |
| edgeR parity | `PYTHON=.venv/bin/python make milo-full-4m-lineage-edger-parity` and `PYTHON=.venv/bin/python make milo-full-4m-lineage-matched-edger-parity` | Exact-neighborhood count-model parity tables. |
| Official MiloR subset sensitivity | `PYTHON=.venv/bin/python make milor-lineage-subset-parity` and `PYTHON=.venv/bin/python make milor-lineage-matched-subset-parity` | Canonical MiloR subset DA summaries. |
| Claim gates | `PYTHON=.venv/bin/python make scvi-scaled-comparison` and `PYTHON=.venv/bin/python make scvi-embedding-claim-gates` | Publication-facing scVI role and marker-concordance gates. |

## Profile C: Final Submission Build

| Check | Command | Required result |
| --- | --- | --- |
| Python checks | `.venv/bin/python -m ruff check .` | Pass. |
| Tests | `PYTHON=.venv/bin/python make test` | Pass. |
| Provenance | `PYTHON=.venv/bin/python make output-provenance` | Zero missing non-deferred outputs. |
| Tables and figures | `PYTHON=.venv/bin/python make publication-tables` and `PYTHON=.venv/bin/python make manuscript-figures` | All manuscript-facing artifacts resolve. |
| Manuscript PDF | `make manuscript` | Pass in an environment with `latexmk` or `pdflatex` plus `bibtex`. |

## Acceptance Notes

- Do not treat Profile A as a full-from-raw reproducibility claim; it refreshes manuscript artifacts from existing inputs.
- Profile B is required before the final reproducibility gate can be marked complete unless the 4M substrate is restored from remote storage with checksums.
- `make manuscript` remains blocked on machines without TeX tooling. That is an environment issue, not an analysis-code failure.

