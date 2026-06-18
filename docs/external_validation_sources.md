# External Validation Source Notes

Updated: 2026-06-15

This note tracks actionable external validation sources for the ORA project. It separates datasets that can support donor-level ORA validation from datasets that are useful only for marker or tissue sanity checks.

## GSE184117: Human Olfactory Aging / Presbyosmia

- Title: Aging-related olfactory loss is associated with olfactory stem cell transcriptional alterations in humans.
- Repository: GEO `GSE184117`; BioProject `PRJNA763212`; SRA `SRP337058`.
- Source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184117
- Supplement: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184117/suppl/GSE184117_RAW.tar
- Series matrix: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184117/matrix/GSE184117_series_matrix.txt.gz
- Cohort: 3 normosmic adults and 3 presbyosmic adults older than 65, plus one culture sample.
- Primary use: external aging/presbyosmia validation for Gateway ORA composition and module signals.
- Immediate adapter work:
  - Download and inspect the raw TAR. Completed locally on 2026-06-15; the archive contains 21 files: matrix, barcodes, and features triplets for 3 controls, 3 presbyosmic samples, and 1 culture sample.
  - Resolve donor/sample metadata, age, and normosmia/presbyosmia status. Completed locally on 2026-06-15 from the GEO series matrix: usable biopsy samples are 3 normosmic controls ages 71, 66, and 51, and 3 presbyosmic/hyposmic/anosmic donors ages 78, 74, and 73; one culture sample from the 51-year-old donor is excluded from biopsy validation.
  - Score sample-level raw 10x modules. Completed locally on 2026-06-15 for 14 curated/published modules across 7 samples; contrasts are descriptive only because n=3 versus n=3 and no public cell labels are present.
  - Score marker-only coarse composition panels. Completed locally on 2026-06-15 for 11 curated olfactory/epithelial/immune marker panels plus an unassigned bin; contrasts are descriptive only and do not replace reference-mapped cell labels.
  - Resolve cell labels or run a defensible reference-mapping/cell-annotation adapter before attempting Gateway-compatible donor composition features.
  - Harmonize any recovered cell labels to Gateway aliases and quantify missing features.

```bash
mkdir -p data/external
curl -L -o data/external/GSE184117_RAW.tar \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184117/suppl/GSE184117_RAW.tar
curl -L -o data/external/GSE184117_series_matrix.txt.gz \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184117/matrix/GSE184117_series_matrix.txt.gz
python scripts/inspect_external_raw_archive.py \
  --archive data/external/GSE184117_RAW.tar \
  --dataset-id oliva_2022
python scripts/score_external_10x_modules.py \
  --archive data/external/GSE184117_RAW.tar \
  --metadata data/external/GSE184117_series_matrix.txt.gz \
  --dataset-id oliva_2022
python scripts/score_external_10x_markers.py \
  --archive data/external/GSE184117_RAW.tar \
  --metadata data/external/GSE184117_series_matrix.txt.gz \
  --dataset-id oliva_2022
python scripts/summarize_external_evidence.py \
  --external-config configs/external_datasets.yaml
```

Current blocker: public GEO metadata resolves age and olfaction status, and raw matrices now support marker-only sample sanity checks, but the supplementary TAR does not include cell labels. Donor-level ORA feature validation still requires author annotations, reference mapping, or a documented de novo annotation workflow.

## GSE151973: Human Olfactory Neuroepithelium Bulk Reference

- Title: A gateway for SARS-CoV-2 infection in the human olfactory neuroepithelium.
- Repository: GEO `GSE151973`; BioProject `PRJNA637909`; SRA `SRP266345`.
- Source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151973
- Primary use: olfactory versus respiratory epithelium marker sanity check.
- Limitation: bulk tissue RNA-seq, so it cannot validate donor-level ORA aging predictions directly.

## Current Priority

1. Implement the `GSE184117` cell-label/reference-mapping adapter first.
2. Use `GSE151973` only as a marker/tissue sanity check unless richer metadata become useful.
3. Continue searching for independent AD/PD olfactory single-cell or spatial datasets; no immediately downloadable donor-level AD/PD olfactory single-cell dataset is confirmed yet.
