# Data And Code Availability Draft

Updated: 2026-06-23

Purpose: hold the submission-ready availability language and the remaining fields that must be replaced at freeze.

## Draft Manuscript Text

The primary source dataset is the Gateway human olfactory epithelial atlas, available through DOI `10.64898/2026.06.10.731272` and the configured CELLxGENE collection ID `8b35aa1f-6bcf-4a51-abc3-a3f336a44ae6`. External public datasets used or assessed for validation include `GSE184117` and `GSE151973`; additional public-search context is documented in the repository validation log.

Analysis code, configuration files, manuscript source, and reproducibility documentation are available at `https://github.com/yashvirsabharwal/olfactory-regenerative-age`. The final submission should cite a frozen release tag or commit SHA: `[INSERT_FINAL_RELEASE_OR_COMMIT]`.

Large raw, processed, model, and result artifacts are intentionally not stored in Git. The repository records command provenance in `configs/command_manifest.yaml` and output provenance in `results/reports/output_provenance.tsv`. Heavyweight manuscript-supporting artifacts, including the full 4M reduced scVI H5AD and associated model weights, are currently identified with SHA-256 checksums in `docs/large_artifact_manifest.md` and staged on `mia.ninds.nih.gov`; before submission, these artifacts should be deposited in durable controlled storage or a public/institutional archive where allowed. The final archival location should be cited here: `[INSERT_ARTIFACT_ARCHIVE_URI]`.

## Freeze Fields

| Field | Current value | Required before submission |
| --- | --- | --- |
| Code repository | `https://github.com/yashvirsabharwal/olfactory-regenerative-age` | Confirm public visibility or reviewer access. |
| Code version | Current moving `main` branch | Replace with release tag or commit SHA. |
| Large artifact archive | Current compute staging on `mia` only | Replace with stable archive URI. |
| Gateway source | DOI `10.64898/2026.06.10.731272`; CELLxGENE ID `8b35aa1f-6bcf-4a51-abc3-a3f336a44ae6` | Confirm final citation and access URL. |
| External data | `GSE184117`, `GSE151973`; context search also logged for `GSE151346`, `GSE290883`, `GSE290884`, `GSE235330`, `GSE324335` | Cite only datasets used in manuscript claims; keep context-only datasets in validation log if not cited. |
| Checksums | `docs/large_artifact_manifest.md` | Refresh after final archive copy. |
| License | TBD | Add repository/software license before submission if absent. |

## Acceptance Gate

Data/code availability is drafted but not submission-complete until the final code version and artifact archive URI are inserted.

