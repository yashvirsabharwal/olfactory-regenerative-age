# External Validation Source Notes

Updated: 2026-06-15

This note tracks actionable external validation sources for the ORA project. It separates datasets that can support donor-level ORA validation from datasets that are useful only for marker or tissue sanity checks.

## GSE184117: Human Olfactory Aging / Presbyosmia

- Title: Aging-related olfactory loss is associated with olfactory stem cell transcriptional alterations in humans.
- Repository: GEO `GSE184117`; BioProject `PRJNA763212`; SRA `SRP337058`.
- Source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184117
- Supplement: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184117/suppl/GSE184117_RAW.tar
- Cohort: 3 normosmic adults and 3 presbyosmic adults older than 65, plus one culture sample.
- Primary use: external aging/presbyosmia validation for Gateway ORA composition and module signals.
- Immediate adapter work:
  - Download and inspect the raw TAR.
  - Resolve donor/sample metadata, age, normosmia/presbyosmia status, and cell labels.
  - Build an AnnData adapter or donor-feature adapter depending on the file structure.
  - Harmonize cell labels to Gateway aliases and quantify missing features.

```bash
mkdir -p data/external
curl -L -o data/external/GSE184117_RAW.tar \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184117/suppl/GSE184117_RAW.tar
```

## GSE151973: Human Olfactory Neuroepithelium Bulk Reference

- Title: A gateway for SARS-CoV-2 infection in the human olfactory neuroepithelium.
- Repository: GEO `GSE151973`; BioProject `PRJNA637909`; SRA `SRP266345`.
- Source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151973
- Primary use: olfactory versus respiratory epithelium marker sanity check.
- Limitation: bulk tissue RNA-seq, so it cannot validate donor-level ORA aging predictions directly.

## Current Priority

1. Implement the `GSE184117` raw adapter first.
2. Use `GSE151973` only as a marker/tissue sanity check unless richer metadata become useful.
3. Continue searching for independent AD/PD olfactory single-cell or spatial datasets; no immediately downloadable donor-level AD/PD olfactory single-cell dataset is confirmed yet.
