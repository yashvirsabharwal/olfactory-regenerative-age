# DESeq2 Parity Decision

Updated: 2026-06-15

## Decision

DESeq2 remains installed in the local R environment, but it is not promoted as a required third genome-wide disease-DE engine for the current manuscript path.

## Rationale

- The disease contrasts are donor-limited: 5 AD and 5 PD donors, with many fine cell states dropping below interpretable case/control balance.
- The strongest DE concern is confounding, not lack of a third count model. Matched FLEX v2/device reruns and sentinel-gene audits directly address that concern.
- edgeR quasi-likelihood and limma-voom now provide two independent pseudobulk engines. limma-voom is substantially more conservative and removes the matched PD sex-linked sentinel signal seen under edgeR.
- Running thousands of DESeq2 fits across fine cell states would add convergence, dispersion, and shrinkage decisions without resolving the main donor-balance and metadata limitations.

## Revisit Criteria

Run DESeq2 parity if any of the following become true:

- Independent AD/PD olfactory validation data provide larger disease/control donor counts.
- Reviewers request DESeq2 specifically for a shortlisted subset of robust edgeR/limma hits.
- The analysis narrows to a small number of predeclared cell states and genes where DESeq2 diagnostics can be reviewed manually.

## Current Reporting Rule

Genome-wide DE biological emphasis should require donor-balance review, sentinel-gene audit, matched FLEX v2/device sensitivity, and edgeR/limma agreement. DESeq2 is documented as deferred rather than absent.
