# Reviewer Risk Memo

Updated: 2026-06-23

Purpose: prepare the paper for likely reviewer objections before submission. This memo is internal and should inform the cover letter, limitations, and response-to-review strategy.

| Risk | Severity | Reviewer concern | Planned response | Supporting artifacts | Remaining action |
| --- | --- | --- | --- | --- | --- |
| External validation is small and mixed. | high | `GSE184117` has only 3 vs 3 biopsy samples and concordance is not uniformly positive. | State that external evidence is small-n support/context, not independent replication; document final search exhaustion and author-label request. | `docs/external_validation_final_search.md`, `docs/gse184117_label_request.md`, `results/tables/manuscript_table_external_validation_strength.tsv`. | Send label request and log outcome. |
| ORA is not an accurate age clock. | high | MAE is modest and calibration is under-dispersed. | Frame ORA as a reproducible, interpretable tissue-state axis above shuffled-age nulls, not an absolute chronological-age clock. | `results/tables/ora_model_card.tsv`, `results/tables/ora_permutation_empirical.tsv`, Figure 3. | Keep title/abstract/results language conservative. |
| Disease projection is confounded. | high | AD/PD groups have 5 donors each and share technical context. | Keep AD/PD projection exploratory; emphasize frozen healthy-trained projection, matched-context diagnostics, and label permutation. | `results/tables/manuscript_table_ndd_guardrails.tsv`, `results/tables/ndd_label_permutation.tsv`. | Avoid disease-biomarker wording. |
| Genome-wide DE is method-sensitive. | medium | edgeR and limma-voom differ; sentinel categories appear in disease contrasts. | Treat DE as audited hypothesis generation; report matched analyses and sentinel audits. | `results/tables/manuscript_table_de_audit_summary.tsv`, Figure 6 / ED Figure 4. | Do not make disease-mechanism claims from genome-wide DE alone. |
| Milo-style analysis is not official MiloR. | medium | Full-scale workflow is Python Milo-style, not canonical MiloR. | Label it Milo-style; report exact-neighborhood edgeR parity and official MiloR subset sensitivity. | `results/tables/milo_full_4m_lineage_edger_parity_summary.tsv`, `results/tables/milor_lineage_subset_summary.tsv`, `docs/run_hierarchy.md`. | Keep Early iOSN as guarded exact-neighborhood subclaim. |
| Large artifacts limit rerunnability. | medium | Full 4M H5AD/model artifacts are too large for Git. | Provide manifest, checksums, rerun profile, remote compute notes, and stable archive URI before submission. | `docs/large_artifact_manifest.md`, `docs/manuscript_rerun_profile.md`, `results/reports/output_provenance.tsv`. | Move `mia` artifacts to durable archive and update URI. |
| Cross-sectional data cannot prove regeneration dynamics. | medium | Composition and neighborhoods may reflect sampling or state, not lineage flux. | Use associational wording; explicitly prohibit lineage-flux claims. | `docs/claim_ledger.md`, `docs/limitations.md`, manuscript limitations. | Final claim-language audit after PDF compile. |
| Manuscript package may fail technical checks. | medium | Missing figures, citations, supplements, or TeX build can delay review. | Use `make manuscript-check`; package currently has 0 asset failures and one local TeX blocker. | `docs/manuscript_package_check.md`, `docs/supplement_inventory.md`. | Compile PDF in TeX-enabled environment. |

## Submission Positioning

Lead with the conservative resource/reanalysis value: a reproducible, donor-level, interpretable healthy olfactory epithelial aging axis. Do not lead with disease projection, external validation strength, or latent-neighborhood mechanism.

