# Supplement Inventory

Updated: 2026-06-23

Purpose: enumerate manuscript-adjacent figures, tables, and documentation that should be bundled or cross-referenced for journal review.

## Main Figures

| Figure | File | Claim gate | Status |
| --- | --- | --- | --- |
| Figure 1 | `results/figures/manuscript_figure1_design.pdf` | Study design and claim gates | Present; referenced in `main.tex`. |
| Figure 2 | `results/figures/manuscript_figure2_age_composition.pdf` | Core ORA composition claim | Present; referenced in `main.tex`. |
| Figure 3 | `results/figures/manuscript_figure3_modeling.pdf` | Core ORA model/null/calibration claim | Present; referenced in `main.tex`. |
| Figure 4 | `results/figures/manuscript_figure4_feature_biology.pdf` | Feature interpretation and module annotation | Present; referenced in `main.tex`. |
| Figure 5 | `results/figures/manuscript_figure5_external_ndd.pdf` | External and NDD guarded extensions | Present; referenced in `main.tex`. |
| Figure 6 | `results/figures/manuscript_figure6_de_latent.pdf` | DE audit and latent-neighborhood guarded extension | Present; referenced in `main.tex`. |

## Extended-Data Figures

| Extended figure | File | Purpose | Status |
| --- | --- | --- | --- |
| ED Figure 1 | `results/figures/extended_data_figure1_model_card.pdf` | ORA model-card detail. | Present. |
| ED Figure 2 | `results/figures/extended_data_figure2_external_evidence.pdf` | External validation evidence ledger. | Present. |
| ED Figure 3 | `results/figures/extended_data_figure3_scvi_validation.pdf` | scVI validation and marker-continuity gates. | Present. |
| ED Figure 4 | `results/figures/extended_data_figure4_de_audit.pdf` | Genome-wide DE audit summaries. | Present. |
| ED Figure 5 | `results/figures/extended_data_figure5_latent_robustness.pdf` | Latent neighborhood robustness and parity. | Present. |
| ED Figure 6 | `results/figures/extended_data_figure6_ndd_guardrails.pdf` | NDD projection guardrails. | Present. |

## Manuscript Tables

| Table | File | Rows | Claim gate | Status |
| --- | --- | ---: | --- | --- |
| Cohort | `results/tables/manuscript_table_cohort.tsv` | 4 | Cohort/design | Present. |
| Model card | `results/tables/manuscript_table_model_card.tsv` | 12 | Core ORA claim | Present. |
| External validation strength | `results/tables/manuscript_table_external_validation_strength.tsv` | 13 | External validation guardrail | Present. |
| Latent neighborhood gates | `results/tables/manuscript_table_latent_neighborhood_gates.tsv` | 10 | Latent/neighborhood support | Present. |
| DE audit summary | `results/tables/manuscript_table_de_audit_summary.tsv` | 8 | Disease/DE guardrail | Present. |
| NDD guardrails | `results/tables/manuscript_table_ndd_guardrails.tsv` | 26 | NDD projection guardrail | Present. |

## Reviewer-Facing Documentation

| Artifact | Role | Submission use |
| --- | --- | --- |
| `docs/journal_acceptance_tracker.md` | Gate dashboard and milestone tracker. | Internal command center; can inform response-to-review strategy. |
| `docs/claim_ledger.md` | Supported, exploratory, deferred, and prohibited claims. | Final claim audit before submission. |
| `docs/external_validation_final_search.md` | Search exhaustion log. | Evidence for maximum validation effort. |
| `docs/gse184117_label_request.md` | Author-label request template/log. | Attach outcome to response strategy; usually not submitted unless requested. |
| `docs/large_artifact_manifest.md` | Large artifact locations and checksums. | Basis for data/code availability and reproducibility supplement. |
| `docs/manuscript_package_check.md` | Citation, figure, table, and TeX audit. | Internal preflight; can be regenerated with `make manuscript-check`. |
| `docs/data_code_availability.md` | Draft data/code availability statement. | Copy into manuscript/submission system after archive URIs are final. |
| `docs/reviewer_risk_memo.md` | Objection-response memo. | Internal cover-letter and response-to-review preparation. |

## Inventory Gate

Supplement inventory is complete for the current manuscript package, except for target-journal formatting decisions and stable archive URIs for heavyweight artifacts.

