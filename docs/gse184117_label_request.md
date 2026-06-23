# GSE184117 Cell-Label Request Log

Updated: 2026-06-23

Purpose: resolve whether the `GSE184117` authors can provide original cell labels or annotation metadata for the external presbyosmia validation. This is a documentation artifact until the request is actually sent.

## Request Status

| Field | Value |
| --- | --- |
| Dataset | `GSE184117` / Oliva et al. |
| Public contact from GEO | Bradley Goldstein, Duke University, `bradley.goldstein@duke.edu` |
| Current request state | Draft ready; not yet sent from project email. |
| Tracker gate | External Validation / M2 |
| Stop condition | Labels integrated; authors decline; files unavailable; or no response after a documented waiting period. |
| Claim rule | Do not promote `GSE184117` beyond small-n supportive/context evidence unless labels or defensible mapping support donor-level ORA feature replication. |

## Requested Materials

| Requested item | Why it matters | Minimum acceptable substitute |
| --- | --- | --- |
| Per-cell author cell-type labels | Tests ORA cell-state features without relying only on marker or reference-transfer labels. | A Seurat/AnnData object with cell annotations. |
| Barcode-to-sample mapping after QC | Aligns annotations to the public raw 10x matrices. | Metadata table with sample prefixes and filtered cell barcodes. |
| Donor/sample phenotype metadata | Confirms age and normosmic/presbyosmic grouping used for contrasts. | Confirmation that GEO sample-level metadata are sufficient. |
| Annotation/QC notes | Documents whether labels are final publication labels and what cells were excluded. | Brief method notes or citation to the exact published annotation procedure. |
| Redistribution/citation guidance | Clarifies whether derived validation tables can be shared in the repository/manuscript supplement. | Permission to cite and summarize without redistributing raw annotation files. |

## Draft Email

Subject: Request for GSE184117 cell annotations for olfactory-aging validation

Dear Dr. Goldstein and colleagues,

I am reanalyzing public human olfactory epithelium single-cell datasets as external validation for a Gateway-atlas study of a healthy-donor olfactory regenerative aging axis. `GSE184117` is the closest public human olfactory aging/presbyosmia dataset we have found. We have used the public GEO raw 10x matrices and series metadata, but I do not see per-cell author annotation tables among the public supplementary files.

Would you be willing to share the original per-cell annotations or metadata for `GSE184117`, ideally including filtered cell barcodes, sample IDs, cell-type labels, QC/filtering notes, and any sample-level age or normosmia/presbyosmia fields that can be cited?

We would use the labels only for conservative external validation/context, cite the GEO record and publication, and avoid diagnostic or age-clock claims. If labels are not available for redistribution, even confirmation of the annotation fields or a non-redistributable metadata file for internal validation would help us keep the manuscript appropriately qualified.

Thank you for considering this request.

Sincerely,

[Sender name]

## Follow-Up Ledger

| Date | Action | Outcome | Tracker update |
| --- | --- | --- | --- |
| 2026-06-23 | Draft request note created. | Ready for user/project email send. | M2 label-recovery item moves from `not started` to `draft ready`. |

