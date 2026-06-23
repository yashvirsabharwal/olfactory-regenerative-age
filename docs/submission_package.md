# Submission Package

Updated: 2026-06-23

Purpose: one journal-submission workspace for administrative fields, data/code availability, cover-letter text, supplement inventory, reviewer risks, and final preflight actions. This replaces the scattered submission admin notes.

## Current Verdict

| Area | Status | Evidence | Remaining action |
| --- | --- | --- | --- |
| Manuscript source | ready except PDF build | `manuscript/main.tex`, `manuscript/references.bib`, `results/reports/manuscript_package_check.md`. | Compile in a TeX-enabled environment. |
| Figures and tables | ready except visual PDF QA | Six main figures, six extended-data figures, six manuscript tables, `results/reports/publication_tables.md`. | Inspect final compiled PDF. |
| Claims | ready with conservative framing | `docs/claim_ledger.md`, `docs/journal_acceptance_tracker.md`. | Final title/abstract/caption audit. |
| External validation | exhausted but label request pending | `docs/external_validation_final_search.md`, `docs/gse184117_label_request.md`, `results/tables/external_candidate_matrix.tsv`. | Send/log GSE184117 label request outcome. |
| Reproducibility | partial | Makefile, command manifest, output provenance, `docs/large_artifact_manifest.md`, `docs/manuscript_rerun_profile.md`. | Add durable artifact archive URI. |
| Author/admin fields | pending | This document. | Fill author-confirmed statements. |

## Final Submission Sequence

| Step | Artifact or command | Pass condition |
| --- | --- | --- |
| 1 | `PYTHON=.venv/bin/python make publication-tables` | Manuscript tables regenerate. |
| 2 | `PYTHON=.venv/bin/python make manuscript-figures` | Main and extended-data figure PDFs regenerate. |
| 3 | `PYTHON=.venv/bin/python make manuscript-check` | Citation and asset checks pass; TeX status is explicit. |
| 4 | `PYTHON=.venv/bin/python make output-provenance` | No missing non-deferred outputs. |
| 5 | `make manuscript` | `manuscript/main.pdf` compiles. |
| 6 | Archive artifacts | Stable URIs and SHA-256 checksums are in `docs/large_artifact_manifest.md`. |
| 7 | Freeze code | Release tag or final commit SHA is recorded below. |
| 8 | Submit | Cover letter, PDF/source, figures, tables, supplement, and availability text are target-journal formatted. |

## Data And Code Availability Draft

The primary source dataset is the Gateway human olfactory epithelial atlas, available through DOI `10.64898/2026.06.10.731272` and the configured CELLxGENE collection ID `8b35aa1f-6bcf-4a51-abc3-a3f336a44ae6`. External public datasets used or assessed for validation include `GSE184117` and `GSE151973`; context-only sources are recorded in `docs/external_validation_final_search.md` and `results/tables/external_candidate_matrix.tsv`.

Analysis code, configuration files, manuscript source, and reproducibility documentation are available at `https://github.com/yashvirsabharwal/olfactory-regenerative-age`. Replace the moving branch reference with a frozen release tag or commit SHA before submission: `[INSERT_FINAL_RELEASE_OR_COMMIT]`.

Large raw, processed, model, and result artifacts are intentionally not stored in Git. The repository records command provenance in `configs/command_manifest.yaml` and output provenance in `results/reports/output_provenance.tsv`. Heavyweight manuscript-supporting artifacts are identified with SHA-256 checksums in `docs/large_artifact_manifest.md`; final submission needs a durable archive URI: `[INSERT_ARTIFACT_ARCHIVE_URI]`.

## Administrative Fields

| Field | Status | Required content |
| --- | --- | --- |
| Target journal and article type | pending | Exact journal, editor if known, resource/reanalysis article category. |
| Final author list | pending | Ordered authors, affiliations, ORCID IDs if required. |
| Corresponding author | pending | Name, email, mailing address if required. |
| Author contributions | pending | CRediT roles or journal-specific format. |
| Funding | pending | Grant numbers, institutional support, or no-specific-funding statement. |
| Competing interests | pending | Author-confirmed conflict statement. |
| Ethics/consent | pending | Source-study human subjects/consent statement or citation-based reuse language. |
| Acknowledgements | pending | Data providers, compute resources, colleagues, non-author contributions. |
| Suggested/opposed reviewers | optional | Names, institutions, emails, expertise, and conflict screening. |
| Preprint decision | pending | Whether and where to post before submission. |
| License | pending | Repository/software license confirmation. |

## Supplement Inventory

| Item | Files | Status |
| --- | --- | --- |
| Main figures | `results/figures/manuscript_figure*.pdf` | Present; final PDF visual QA needed. |
| Extended-data figures | `results/figures/extended_data_figure*.pdf` | Present; target journal may require supplement-only formatting. |
| Manuscript tables | `results/tables/manuscript_table_*.tsv` | Present and indexed in `results/reports/publication_tables.md`. |
| Supplement tables | Pipeline result TSVs listed through `results/reports/publication_tables.md` and provenance. | Select target-journal bundle after venue choice. |
| Reproducibility supplement | `docs/run_hierarchy.md`, `docs/manuscript_rerun_profile.md`, `docs/large_artifact_manifest.md`. | Needs archive URI. |

## Reviewer Risk Memo

| Risk | Severity | Planned response | Supporting artifacts | Remaining action |
| --- | --- | --- | --- | --- |
| External validation is small and mixed. | high | Treat GSE184117 as small-n support/context, not independent replication. | Final search log, label request, validation-strength table. | Send/log label request outcome. |
| ORA is not an accurate age clock. | high | Frame ORA as a reproducible tissue-state axis above shuffled-age nulls. | Model card, permutation null, calibration, Figure 3. | Keep title/abstract conservative. |
| Disease projection is confounded. | high | Keep AD/PD projection exploratory and guarded by matched-context diagnostics. | NDD guardrail table, label permutation. | Avoid biomarker wording. |
| Genome-wide DE is method-sensitive. | medium | Treat DE as audited hypothesis generation. | DE audit table, edgeR/limma/matched/sentinel outputs. | Do not make disease-mechanism claims from DE alone. |
| Milo-style analysis is not official MiloR. | medium | Label full-scale analysis Milo-style; use edgeR parity and official MiloR subset sensitivity. | Milo parity tables, `docs/run_hierarchy.md`. | Keep Early iOSN guarded. |
| Large artifacts limit rerunnability. | medium | Provide manifest, checksums, rerun profile, and stable archive URI. | Artifact manifest, provenance, rerun profile. | Archive final H5AD/model artifacts. |
| Cross-sectional data cannot prove regeneration dynamics. | medium | Use associational wording and prohibit lineage-flux claims. | Claim ledger, SOTA tracker, manuscript limitations. | Final claim-language audit. |

## Cover Letter Skeleton

Dear `[Editor name]`,

We are pleased to submit our manuscript, "An Interpretable Olfactory Regenerative Aging Axis from Human Single-Cell Epithelial Composition," for consideration as a computational single-cell resource article in `[Journal]`.

This study reanalyzes the Gateway 4M-cell human olfactory epithelial atlas to define a donor-level olfactory regenerative aging axis in healthy donors. The manuscript emphasizes a conservative and reproducible claim: healthy human olfactory epithelial composition contains a modest but interpretable age-associated regenerative-state signal. We do not frame ORA as an absolute biological-age clock or disease diagnostic model.

The work combines donor-level composition and module features, repeated donor-level cross-validation, shuffled-age null testing, calibration diagnostics, feature interpretation, external presbyosmia-context validation, AD/PD projection guardrails, genome-wide DE audits, and full-scale scVI/Milo-style latent-neighborhood analyses with edgeR and official MiloR sensitivity checks.

The manuscript is original, is not under consideration elsewhere, and all authors have approved the submission. `[INSERT CONFLICT OF INTEREST STATEMENT]`

Sincerely,

`[Corresponding author name and contact]`
