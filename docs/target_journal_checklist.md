# Target Journal Checklist

Updated: 2026-06-23

Target strategy: specialty computational biology / single-cell resource journal. This document should be converted to a journal-specific checklist after the target venue is chosen.

## Generic Specialty-Journal Requirements

| Requirement | Current status | Evidence | Remaining action |
| --- | --- | --- | --- |
| Clear article category | ready | Framed as computational single-cell resource/reanalysis. | Pick exact journal article type. |
| Conservative main claim | ready | `docs/claim_ledger.md`, `docs/journal_acceptance_tracker.md`. | Final abstract/title audit. |
| Reproducible code | partial | GitHub repo, Makefile, command manifest, tests. | Freeze release/tag and confirm license. |
| Source data identifiers | ready | Gateway DOI/CELLxGENE ID, GEO accessions in availability draft. | Confirm final citation formatting. |
| Large artifact policy | partial | `docs/large_artifact_manifest.md` with local and `mia` checksums. | Add stable archival URI. |
| External validation effort | mostly ready | Final search log and request draft. | Send/log GSE184117 label request outcome. |
| Statistical guardrails | ready | Repeated CV, permutation, calibration, DE audits, NDD guardrails. | Keep claims aligned in final text. |
| Figures and tables | ready except PDF QA | `docs/manuscript_package_check.md`, `docs/supplement_inventory.md`. | Compile and inspect final PDF. |
| Supplement inventory | ready | `docs/supplement_inventory.md`. | Format to target journal. |
| Reporting checklist | pending | No target-specific checklist yet. | Complete after journal selection. |
| Ethics/human subjects language | pending | Gateway source citation exists. | Add source-study ethics/consent statement or citation-based language. |
| Author contributions | pending | Not drafted. | Add CRediT roles. |
| Funding/competing interests | pending | Placeholder only. | Add author-confirmed statements. |
| Administrative fields | pending | `docs/author_submission_fields.md`. | Finalize author-confirmed metadata in the journal system. |

## Decision Gate

The project can move from "publication-ready internally" to "submission-ready externally" only after:

- A target journal is selected.
- TeX PDF compile passes.
- Stable artifact archive URI is inserted.
- GSE184117 label request is sent/logged or explicitly deferred with limitation.
- Author/funding/ethics/conflict statements are finalized.
