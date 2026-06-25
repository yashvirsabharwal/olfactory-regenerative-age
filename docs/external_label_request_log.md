# External Label Request Log

Updated: 2026-06-25

## GSE184117 Author Label Request

Status: Response received; no separate per-cell manual annotation table available.  
Dataset: GSE184117, "Aging-related olfactory loss is associated with olfactory stem cell transcriptional alterations in humans."  
Public source checked: NCBI GEO accession page, accessed 2026-06-24.  
GEO contact: Bradley Goldstein, Duke University, `bradley.goldstein@duke.edu`.  
Sent date: 2026-06-25, reported by project owner.  
Response date: 2026-06-25, reported by project owner.  
Response summary: Raw sequencing files and 10x matrix files are available through GEO. Sample names and clinical metadata in manuscript Table 1 can be matched to the corresponding GEO record name. The authors did not deposit a separate per-cell manual annotation table, so those labels are not available as a separate file. Public matrices allow independent reanalysis and annotation, and donor/sample metadata can be matched using GEO records and manuscript tables.

### Materials Requested

- Per-cell annotation labels used in the publication or internal analysis.
- Barcode-to-sample/GSM mapping for the raw 10x matrices.
- Donor/sample phenotype metadata, including olfactory status, smell-test score, age if shareable, sex if shareable, and biopsy/culture flag.
- QC notes, annotation notes, marker definitions, and exclusion criteria.
- Redistribution and citation guidance for derived labels or metadata.

### Draft Email

To: `bradley.goldstein@duke.edu`  
Subject: Request for GSE184117 olfactory epithelium cell labels and sample metadata

```text
Dear Dr. Goldstein,

I am reanalyzing public human olfactory epithelium single-cell datasets for a manuscript on donor-level regenerative aging signatures. GSE184117 is the most directly relevant public olfactory aging/presbyosmia dataset we have found, and we are using it as a guarded external validation/context dataset.

The GEO record provides the raw 10x-style matrices and sample-level descriptions, but we have not found the per-cell labels and detailed barcode/sample annotation needed to compare donor-level cell-state features cleanly against our Gateway olfactory atlas analysis.

Would you be willing to share any of the following materials, if available?

1. Per-cell annotation labels used in the publication or internal analysis, ideally with barcode, sample/GSM, cell type or state label, and QC inclusion status.
2. Barcode-to-sample/GSM mapping for the raw matrices.
3. Donor/sample phenotype metadata, including olfactory status, smell-test score, age if shareable, sex if shareable, and biopsy/culture flag.
4. QC notes, annotation notes, marker definitions, and exclusion criteria.
5. Guidance on redistribution and citation, including whether derived labels/metadata may be included in a supplementary table or should only be used internally for validation.

We would cite the original study and GEO accession prominently and follow any redistribution limits you specify. If these labels are unavailable or cannot be shared, a brief confirmation would still be very helpful so we can accurately document the external-validation limitation.

Thank you for considering this request.

Best,
[sender name]
[institution / project]
[reply email]
```

### Outcome Log

| Date | Event | Recipient | Sender | Outcome | Follow-up |
| --- | --- | --- | --- | --- | --- |
| 2026-06-24 | Draft prepared from GEO contact record. | `bradley.goldstein@duke.edu` | pending at draft time | Superseded by sent request on 2026-06-25. | None. |
| 2026-06-25 | Label request sent. | `bradley.goldstein@duke.edu` | project owner | Sent; response received the same day. | Superseded by response row below. |
| 2026-06-25 | Author response received. | project owner | Bradley Goldstein / study team | No separate per-cell manual annotation table is available; raw sequencing files and 10x matrices are public through GEO; sample metadata can be matched from GEO record names and manuscript Table 1. | Treat GSE184117 as public-matrix reanalysis with sample-level metadata matching, not author-label replication. |

### Evidence Classification After Response

- Per-cell author labels: unavailable as a separate deposited/shareable file.
- Barcode/sample mapping: use public 10x matrix file organization and GEO record names.
- Donor/sample metadata: match GEO record names to manuscript Table 1.
- Validation status: guarded small-n public reanalysis/contextual validation; not an independent author-label replication.
