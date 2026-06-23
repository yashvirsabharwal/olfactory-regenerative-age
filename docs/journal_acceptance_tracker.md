# Journal Acceptance Milestone Tracker

Updated: 2026-06-23

This document is the journal-readiness command center for moving ORA from a strong internal analysis to a specialty computational biology or single-cell resource submission. It converts the SOTA research tracker into acceptance gates: every likely reviewer objection needs a documented response, supporting artifact, claim decision, and pass/fail status.

Acceptance cannot be guaranteed. In this tracker, "publication ready" means the manuscript is maximally defensible for journal review: the core claim is supported, limitations are explicit, validation has been pushed as far as feasible, and reproducibility artifacts are packaged well enough for review.

## Submission Target And Acceptance Bar

Target article type: specialty journal / computational single-cell resource article.

Working frame: a reproducible and interpretable donor-level reanalysis of the Gateway human olfactory epithelial atlas, centered on a healthy-donor olfactory regenerative aging axis.

Required standard for submission:

- Reproducible atlas reanalysis with clear source data identifiers, command manifest, output provenance, and large-artifact policy.
- Defensible statistics with donor-level splitting, shuffled-age null testing, calibration, sensitivity analysis, and explicit covariate/technical guardrails.
- Transparent limitations for external validation, AD/PD sample size, technical confounding, disease DE, and latent-neighborhood interpretation.
- Maximum feasible external-validation effort before submission, with all validation evidence classified by strength.
- Manuscript, figures, tables, supplement inventory, code availability, and data availability frozen into a reviewable package.

Prohibited claims:

- ORA is an absolute biological-age clock.
- ORA has AD/PD diagnostic utility.
- Cross-sectional data measure lineage flux.
- UMAP-derived structure supports trajectory or neighborhood biology.

Source-of-truth documents:

- `docs/sota_research_tracker.md`
- `docs/claim_ledger.md`
- `docs/run_hierarchy.md`
- `docs/publication_tables.md`
- `manuscript/main.tex`
- `results/reports/output_provenance.tsv`

## Claim Promotion Rules

| Claim class | Promotion rule | Allowed manuscript placement | Current decision |
| --- | --- | --- | --- |
| Healthy ORA aging axis | Gate is complete when repeated CV, shuffled-age null, calibration, sensitivity, and feature interpretation are final. | Main text primary claim. | Promote as modest, reproducible, interpretable tissue-state axis. |
| Module augmentation | Gate is complete when predictive gain and overlapping intervals are reported. | Main text or supplement as biological annotation, not primary performance claim. | Keep as modest annotation/support. |
| GSE184117 external validation | Gate is complete only after maximum label recovery/search effort is documented and evidence is classified. | Main text guarded support plus extended data. | Do not promote beyond small-n supportive/context evidence. |
| AD/PD ORA projection | Gate is complete when all projection, matched-context, uncertainty, and label-permutation guardrails are reported. | Main text exploratory section or supplement. | Keep exploratory only. |
| Genome-wide disease DE | Gate is complete when edgeR, limma-voom, matched context, donor-balance, and sentinel audits agree for any emphasized claim. | Supplement or guarded main-text audit. | Keep hypothesis-generating. |
| Full 4M latent neighborhoods | Gate is complete when full 4M scVI, Python Milo-style, edgeR parity, official MiloR subset sensitivity, age-bin robustness, and program scoring are summarized consistently. | Secondary mechanistic layer. | Promote broad lineage-neighborhood remodeling only; keep Early iOSN as guarded exact-neighborhood subclaim. |
| Pseudotime, cNMF, ligand-receptor, spatial validation | Gate requires new method-specific workflows and sensitivity checks. | Future work unless completed. | Defer. |

## Acceptance Gate Dashboard

| Gate | Status | Owner role | Acceptance criterion | Current evidence and artifacts | Reviewer risk | Claim unlocked | Next action |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Core Claim | complete | Analysis lead | Repeated CV, shuffled-age null, calibration, feature interpretation, and sensitivity outputs are final and reflected in the manuscript. | `results/tables/ora_repeated_cv_summary.tsv`, `results/tables/ora_permutation_empirical.tsv`, `results/tables/ora_calibration.tsv`, `results/tables/ora_feature_interpretation.tsv`, `results/tables/ora_sensitivity_performance.tsv`, `manuscript/main.tex`. | Reviewer may argue the age signal is modest. | Healthy-donor ORA is a modest, reproducible tissue-state axis. | Keep language focused on interpretability and null-model separation, not age-clock accuracy. |
| External Validation | needs work | Validation lead | GSE184117 labels are requested/recovered or documented unavailable; GSE151973 context is assessed; one final public-data search is logged; all evidence is classified as direct, mapped, marker-only, context-only, or blocked. | Current evidence: `docs/external_validation_final_search.md`, `docs/gse184117_label_request.md`, `results/tables/external_validation_evidence.tsv`, `results/tables/external_scanvi_feature_concordance.tsv`, `results/tables/manuscript_table_external_validation_strength.tsv`. | Largest scientific weakness; small n and mixed concordance can limit acceptance. | At most guarded external support unless stronger data are recovered. | Send GSE184117 label request from project email, then update request outcome and manuscript limitations. |
| Latent/Neighborhood Support | complete | Latent analysis lead | Full 4M scVI/Milo-style outputs, edgeR parity, official MiloR subset sensitivity, age-bin robustness, and program enrichment are summarized with conservative claim language. | `results/tables/scvi_embedding_claim_gates.tsv`, `results/tables/milo_full_4m_lineage_summary.tsv`, `results/tables/milo_full_4m_lineage_matched_summary.tsv`, `results/tables/milo_full_4m_lineage_edger_parity_summary.tsv`, `results/tables/milor_lineage_subset_summary.tsv`, `results/tables/milo_full_4m_lineage_matched_program_summary.tsv`, `docs/run_hierarchy.md`. | Reviewer may object that Python Milo-style is not official MiloR or that Early iOSN is narrow. | Secondary broad lineage-neighborhood remodeling. | Keep official MiloR as sensitivity and Early iOSN as guarded exact-neighborhood result. |
| Disease/DE Guardrails | defer with limitation | Disease analysis lead | AD/PD projection remains exploratory; disease DE claims require edgeR/limma/matched/sentinel agreement or stay hypothesis-generating. | `results/tables/ndd_ora_projection_summary.tsv`, `results/tables/ndd_label_permutation.tsv`, `results/tables/manuscript_table_ndd_guardrails.tsv`, `results/tables/manuscript_table_de_audit_summary.tsv`, `docs/claim_ledger.md`. | Disease cohorts are small, technically confounded, and method-sensitive. | No disease-biomarker claim; only exploratory disease observations. | Ensure all disease and DE language stays explicitly exploratory in main text, captions, and supplement. |
| Reproducibility | needs work | Reproducibility lead | Heavyweight artifact locations/checksums are documented; stale workflow paths are removed; manuscript rerun profile exists; provenance reports zero missing non-deferred outputs. | Current non-deferred provenance exists in `results/reports/output_provenance.tsv`; Makefile is the active workflow; `docs/large_artifact_manifest.md` and `docs/manuscript_rerun_profile.md` now document local checksums, `mia` remote checksums for the primary 4M artifacts, and rerun order. | Reviewer may question rerunnability of large latent/neighborhood results. | Engineering claim that outputs are reproducible from documented source state. | Move final 4M H5AD/model artifacts from `mia` to stable archival storage, add stable URIs, and rerun provenance before submission. |
| Manuscript Package | needs work | Manuscript lead | PDF compiles; all figures/tables resolve; citation keys pass; claims match ledger; supplement inventory is complete. | `manuscript/main.tex`, `manuscript/references.bib`, `docs/manuscript_package_check.md`, `docs/submission_package.md`, `results/figures/manuscript_figure*.pdf`, `results/figures/extended_data_figure*.pdf`, `docs/publication_tables.md`. | Local TeX tooling is missing; unresolved figure/table paths would block submission. | Reviewable manuscript package once PDF compiles. | Build in a TeX-enabled environment, visually inspect PDF, and run final claim-language audit. |
| Submission Compliance | needs work | Submission lead | Target journal checklist is completed, cover letter is drafted, data/code availability is finalized, and reviewer-risk memo is written. | `docs/submission_package.md`. | Administrative mismatch can delay or desk-reject an otherwise strong paper. | Submission-ready package after target-specific fields, archive URI, and author statements are finalized. | Pick target journal, add stable artifact URI, finalize author/funding/ethics/conflict fields, and tailor cover letter. |

## Milestone Roadmap

| Milestone | Status | Acceptance criterion | Required artifacts | Claim decision | Dependencies |
| --- | --- | --- | --- | --- | --- |
| M1: Core ORA Claim Lock | complete | Repeated CV, shuffled-age null, calibration, feature interpretation, and sensitivity outputs are final and reflected in manuscript. | ORA model card, repeated CV, permutation, calibration, feature interpretation, sensitivity tables, manuscript Results/Methods text. | Main primary claim is allowed. | None. |
| M2: Maximum External Validation Push | needs work | GSE184117 labels are requested/recovered or documented unavailable; GSE151973 context is assessed; final public-data search is logged; all validation evidence is classified. | `docs/gse184117_label_request.md`, `docs/external_validation_final_search.md`, updated evidence ledger, validation-strength table, manuscript limitation text. | External validation remains guarded unless stronger donor-level evidence is recovered. | Access to public metadata and any author response. |
| M3: Latent And Neighborhood Claim Lock | complete | Full 4M scVI/Milo-style outputs, edgeR parity, official MiloR subset sensitivity, age-bin robustness, and program enrichment are summarized with conservative language. | Run hierarchy, embedding claim gates, Milo summaries, edgeR parity summaries, official MiloR summaries, program summaries, Figure 6 / Extended Data. | Secondary broad lineage-neighborhood claim is allowed; Early iOSN stays guarded. | Large-artifact manifest still needed for reproducibility gate. |
| M4: Disease And DE Guardrail Lock | defer with limitation | AD/PD projection remains exploratory; DE claims require edgeR/limma/matched/sentinel agreement or stay hypothesis-generating. | NDD guardrails table, label permutation, uncertainty table, DE audit table, donor-balance/sentinel summaries, limitations text. | No disease diagnostic or disease-mechanism claim. | None for conservative submission. |
| M5: Reproducibility Package | needs work | Heavyweight artifact locations/checksums are documented; stale workflow paths are removed; manuscript rerun profile exists; provenance reports zero missing non-deferred outputs. | `docs/large_artifact_manifest.md`, `docs/manuscript_rerun_profile.md`, updated command/provenance tables. | Supports reproducibility/resource claim. | Needs stable archival URIs for the full 4M artifacts currently identified on `mia`. |
| M6: Manuscript And Figure Freeze | needs work | PDF compiles; all figures/tables resolve; citation keys pass; claims match ledger; supplement inventory is complete. | Compiled PDF, final figures, final publication tables, `docs/submission_package.md`, `docs/manuscript_package_check.md`, claim-ledger diff check. | Locks text for submission. | TeX environment and final claim-language audit. |
| M7: Submission Readiness | needs work | Target journal checklist is complete, cover letter drafted, data/code availability finalized, and reviewer-risk memo written. | `docs/submission_package.md`, final response strategy. | Ready to submit after target-specific fields are filled. | M2 request outcome, M5 archive URI, M6 PDF compile, author/funding/ethics/conflict statements. |

## Maximum Validation Work Queue

Because the chosen validation stance is maximum validation before submission, M2 has the highest remaining scientific priority.

| Task | Status | Acceptance output | Stop condition |
| --- | --- | --- | --- |
| Request GSE184117 author labels or annotation files | draft ready | `docs/gse184117_label_request.md` or equivalent logged note; tracker updated with request date and outcome. | Labels integrated, authors decline, no response after a documented waiting period, or files are unavailable. |
| Reassess GSE151973 as context-only bulk/deconvolution evidence | complete | Note in `docs/external_validation_final_search.md`; evidence ledger updated if useful. | It either clarifies marker specificity as context-only or is documented as not useful. |
| Run one final public search for human olfactory/nasal aging, presbyosmia, anosmia, COVID, AD/PD, spatial, and histology validation sources | complete | Search log with date, queries, candidates, inclusion/exclusion decisions, and registry updates. | No stronger donor-level validation source is found, or a new candidate is promoted to adapter work. |
| Reclassify all validation evidence after the final search | complete | Updated `results/tables/external_validation_evidence.tsv` and `results/tables/manuscript_table_external_validation_strength.tsv`. | Every row has one of: direct, mapped, marker-only, context-only, blocked. |
| Update manuscript language after validation push | pending | Main text, limitations, and captions match the final evidence strength. | No validation language exceeds `docs/claim_ledger.md`. |

## Reproducibility Package Work Queue

| Task | Status | Acceptance output | Notes |
| --- | --- | --- | --- |
| Large-artifact manifest | remote located | Document path, storage location, size, checksum policy, and regeneration command for full 4M H5AD/model artifacts. | Must cover deferred full scVI and Milo-related inputs not tracked by normal provenance. |
| Rerun profile | drafted | A Make target or documented staged profile that regenerates manuscript-facing artifacts without guessing command order. | Can reference remote compute for full 4M steps. |
| Stale workflow cleanup | complete | Makefile is canonical and the stale Snakemake MVP workflow has been removed. | Avoids stale workflow confusion during review. |
| Provenance refresh | current | `PYTHON=.venv/bin/python make output-provenance` reports zero missing non-deferred outputs after final figures/tables. | Deferred artifacts are covered by `docs/large_artifact_manifest.md`. |
| Environment verification | partial | Python and R package checks logged; TeX and MiloR local paths fixed or documented. | Official MiloR rerun path should not point to a missing executable. |

## Manuscript Freeze Checklist

| Check | Required result | Status |
| --- | --- | --- |
| Claim ledger alignment | No main-text, caption, or abstract claim exceeds `docs/claim_ledger.md`. | needs final pass |
| Figure resolution/path check | All main and extended-data figure PDFs exist and are referenced correctly. | pass in `docs/manuscript_package_check.md` |
| Publication table inventory | `docs/publication_tables.md` lists all manuscript-ready tables. | current |
| Citation check | Every citation key in `manuscript/main.tex` exists in `manuscript/references.bib`. | pass in `docs/manuscript_package_check.md` |
| PDF build | `make manuscript` succeeds in a TeX-enabled environment. | blocked locally by missing TeX |
| Supplement inventory | All supplementary tables/figures are named, described, and claim-gated. | drafted in `docs/submission_package.md` |
| Data/code availability | Source data accessions, repository state, generated artifact policy, and large-artifact access are explicit. | drafted; needs final archive URI and release SHA |

## Reviewer Objection Register

| Likely objection | Planned response | Supporting artifacts | Current status | Claim consequence |
| --- | --- | --- | --- | --- |
| The age model is weak. | ORA is framed as a modest but reproducible tissue-state axis, not an absolute age clock; null testing and repeated CV show signal above shuffled age. | `ora_repeated_cv_summary.tsv`, `ora_permutation_empirical.tsv`, `ora_model_card.tsv`, Figure 3. | response ready | Keep primary claim modest. |
| External validation is underpowered. | Push maximum validation, request/recover labels if possible, classify all evidence by strength, and state that external support is small-n/contextual. | `external_validation_evidence.tsv`, `manuscript_table_external_validation_strength.tsv`, final search log, GSE184117 request note. | needs work | No definitive validation claim unless stronger data emerge. |
| Disease cohorts are confounded. | AD/PD projection is frozen-model, small-n, matched-context guarded, and exploratory only. | `ndd_ora_projection_summary.tsv`, `ndd_label_permutation.tsv`, `manuscript_table_ndd_guardrails.tsv`. | response ready | No disease diagnostic or mechanistic claim. |
| Genome-wide DE is method-sensitive. | Report edgeR/limma parity, matched subset, donor-balance, and sentinel audits; use DE as hypothesis generation only. | `manuscript_table_de_audit_summary.tsv`, DE summary and audit tables. | response ready | No primary disease-DE biology claim. |
| Milo-style is not official MiloR. | Primary full-scale map is labeled Python Milo-style; exact-neighborhood edgeR parity and official MiloR subset sensitivity are reported transparently. | `milo_full_4m_lineage_edger_parity_summary.tsv`, `milor_lineage_subset_summary.tsv`, `docs/run_hierarchy.md`. | response ready | Broad lineage-neighborhood remodeling only; Early iOSN remains guarded. |
| Reproducibility depends on large artifacts. | Large artifacts are identified locally or on `mia` with checksums, but still need stable archival URIs before submission. | `docs/large_artifact_manifest.md`, `output_provenance.tsv`, command manifest, rerun profile. | needs archive URI | Resource/reproducibility claim remains incomplete until packaged. |
| The analysis is cross-sectional. | State that ORA and neighborhoods are associational and cannot measure lineage flux. | `docs/claim_ledger.md`, manuscript limitations. | response ready | No lineage-flux claim. |
| Module scores are not UCell/rank based. | Report average log1p module scoring as a deliberate, transparent summary; avoid claiming exact UCell parity. | `module_gene_coverage.tsv`, `module_score_summary.tsv`, Methods. | response ready | Module evidence is biological annotation, not standalone mechanism. |

## Final Submission Readiness Definition

The project is ready for specialty-journal submission when:

- M1, M3, M4, M5, M6, and M7 are complete, and M2 is either complete or explicitly closed with documented validation exhaustion.
- Every claim in the abstract, main text, captions, and supplement maps to a completed gate or is explicitly marked exploratory/limited.
- The compiled PDF, figures, tables, provenance, command manifest, data/code availability statement, cover letter, and reviewer-risk memo are present.
- No known reviewer objection lacks either a supporting artifact or an explicit limitation.
