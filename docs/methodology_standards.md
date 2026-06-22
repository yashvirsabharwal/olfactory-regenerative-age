# ORA Methodology Standards

Updated: 2026-06-22

This document defines the analysis standard for manuscript-facing results. Smaller runs can be used for engineering checks, but manuscript claims should use the strongest available real-data run and the claim gates below.

## Cohort And Splits

- Train ORA age models only on healthy donors with valid age.
- Keep AD/PD donors fully held out for frozen-model projection only.
- Cross-validation is donor-level, never cell-level.
- Technical variables such as chemistry, collection method, sex, site, and yield are diagnostics, covariates, or sensitivity variables; they are not primary biological ORA features.

## Feature And Model Standards

- Primary ORA features are donor-level cell-state proportions, CLR composition features, lineage ratios, and curated module features.
- Model claims require repeated donor-level CV, shuffled-age null checks, calibration diagnostics, and residual stratification.
- Accuracy should be described as modest and reproducible, not as an absolute biological-age clock.
- Exploratory model families are acceptable only under nested donor-level CV and must be reported as exploratory if they reduce interpretability.

## Latent Atlas Standards

- UMAP is visualization-only and cannot support trajectory, Milo, or cNMF claims.
- The current primary latent substrate is the all-cell reduced Gateway scVI model: 4,028,275 cells, 3,000/3,003 HVG-marker genes, finite 10-dimensional `X_scvi`.
- Stratified 250k and lineage-focused 100k scVI models are sensitivity anchors and engineering checks, not final scientific endpoints.
- Latent claims require label purity, technical-mixing, marker-continuity, and seed/sensitivity checks.

## Milo-Style Neighborhood Standards

- Manuscript-facing neighborhood results should use full 4M scVI targets unless explicitly labeled as sensitivity/pilot.
- Current full 4M standards: 20,000 stratified seed neighborhoods, 100 nearest cells, at least 30 donors, donor-level logit-fraction regression, and adjustment for sex, chemistry, collection method, and donor yield.
- Because neighborhoods overlap, significant-neighborhood counts are maps of recurring local signal, not independent discovery counts.
- Full 4M Milo-style results can support a conservative secondary mechanistic layer. Matched FLEX v2/device sensitivity, ORA-theme annotation, curated lineage-neighborhood program scoring, age-bin robustness, exact-neighborhood edgeR parity, and official MiloR subset sensitivity are now complete and should govern claim language.
- Use "Milo-style" for the Python full-scale workflow. Official MiloR was run as a subset sensitivity with independently constructed neighborhoods; it confirms broad age-associated lineage-neighborhood structure but does not independently reproduce matched Early iOSN as the dominant signal.

## External Validation Standards

- GSE184117 is a small-n validation/sanity dataset, not independent proof of ORA.
- Marker-only, mapped-feature, and scANVI/scArches concordance results must be reported as small-n and claim-gated unless a larger independent donor-level dataset is added.
- GSE151973 remains marker-context or deconvolution context, not donor-level aging validation.

## DE And Disease Standards

- Genome-wide DE must be interpreted through donor-balance, sentinel-gene, matched-subset, and edgeR/limma parity audits.
- AD/PD ORA projection is hypothesis-generating only because disease groups have 5 donors each and share FLEX v2/device context.
- No diagnostic or disease-biomarker claim is allowed without independent disease validation.

## Manuscript Claim Hierarchy

- Primary: healthy olfactory epithelial composition/module features encode a modest, reproducible regenerative aging axis.
- Secondary/mechanistic: full-scale scVI neighborhoods suggest age-associated shifts in lineage neighborhoods. Exact Python-neighborhood analyses support a matched Early iOSN/iOSN depletion signal, while official MiloR subset sensitivity supports broader age-associated HBC/sustentacular/suprabasal structure and narrows the Early iOSN language.
- Exploratory: GSE184117 concordance, AD/PD projection, and genome-wide disease DE.
- Deferred: pseudotime, cNMF, CellRank, ligand-receptor mechanism, spatial/histology validation.
