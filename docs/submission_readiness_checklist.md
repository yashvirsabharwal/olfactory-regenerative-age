# Submission Readiness Checklist

Updated: 2026-06-23

Scope: specialty computational biology / single-cell resource submission. This checklist remains journal-generic until a final target venue is selected.

## Current Verdict

| Area | Status | Evidence | Remaining action |
| --- | --- | --- | --- |
| Manuscript source | pass | `manuscript/main.tex` exists and is checked by `make manuscript-check`. | Compile PDF in TeX-enabled environment. |
| Bibliography | pass | `make manuscript-check` found 15 cited keys, 15 bibliography keys, 0 missing, 0 unused. | Re-run after any citation edits. |
| Main figures | pass | Six `manuscript_figure*.pdf` files exist and are referenced in `main.tex`. | Visual QA after final PDF compile. |
| Extended-data figures | pass | Six `extended_data_figure*.pdf` files exist. | Decide whether target journal allows extended-data figures or wants supplement-only figures. |
| Manuscript tables | pass | Six `manuscript_table_*.tsv` files exist and are indexed in `docs/publication_tables.md`. | Attach as supplement or convert to journal table format. |
| TeX/PDF build | blocked locally | `make manuscript-check` reports no local TeX engine. | Run `make manuscript` on a TeX-enabled machine or install `latexmk`/`pdflatex`+`bibtex`. |
| Supplement inventory | drafted | `docs/supplement_inventory.md`. | Convert to target journal supplement template after venue selection. |
| Data/code availability | drafted | `docs/data_code_availability.md`. | Replace `mia` compute paths with stable archival URIs and final commit SHA. |
| Reviewer-risk memo | drafted | `docs/reviewer_risk_memo.md`. | Update after GSE184117 author-label request outcome. |
| Cover letter | drafted | `docs/cover_letter_draft.md`. | Customize to target journal/editor and final author list. |
| Author/admin fields | pending | `docs/author_submission_fields.md`. | Fill final author, contribution, funding, ethics, conflict, and acknowledgement statements. |
| External validation | exhausted but request pending | `docs/external_validation_final_search.md`, `docs/gse184117_label_request.md`. | Send request and log outcome or close with documented no-response window. |
| Large artifacts | remote located | `docs/large_artifact_manifest.md` records local and `mia` checksums. | Archive final H5AD/model artifacts in durable storage and add URIs. |

## Journal-Specific Fields To Fill

| Field | Current value | Required before submission |
| --- | --- | --- |
| Target journal | TBD specialty computational biology / single-cell resource venue | Exact journal name and article type. |
| Word limit | TBD | Confirm abstract, main text, methods, figure legends, supplement limits. |
| Figure format | PDF currently available; PNG also generated | Confirm accepted file types, DPI, color mode, panel labeling. |
| Data policy | Drafted; large artifacts staged on `mia` | Stable archive link and repository commit. |
| Code policy | GitHub repository available | Final public URL, release tag, and license statement. |
| Reporting checklist | Generic single-cell resource expectations | Target-specific checklist if required. |
| Competing interests | Not drafted | Author-confirmed statement. |
| Funding/acknowledgements | Placeholder in manuscript | Author-confirmed statement. |
| Author contributions | Not drafted | CRediT roles or journal-specific format. |
| Ethics/consent | Not drafted | Source-atlas ethics statement or citation-based statement. |

## Final Submission Sequence

| Step | Command or artifact | Pass condition |
| --- | --- | --- |
| 1 | `PYTHON=.venv/bin/python make publication-tables` | Six manuscript tables regenerated. |
| 2 | `PYTHON=.venv/bin/python make manuscript-figures` | Six main and six extended-data figure PDFs regenerated. |
| 3 | `PYTHON=.venv/bin/python make manuscript-check` | 0 failures; TeX either pass or explicitly blocked on local machine. |
| 4 | `PYTHON=.venv/bin/python make output-provenance` | 0 missing non-deferred outputs. |
| 5 | `make manuscript` | `manuscript/main.pdf` compiles in TeX-enabled environment. |
| 6 | Archive artifacts | Stable URIs and SHA-256 checksums added to `docs/large_artifact_manifest.md`. |
| 7 | Freeze repository | GitHub release/tag or final commit SHA recorded in availability statement. |
| 8 | Submit package | Cover letter, manuscript PDF, source files, figures, tables, supplement inventory, availability statement, reviewer-risk memo ready. |
