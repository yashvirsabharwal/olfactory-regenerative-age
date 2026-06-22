# Milo-Style Implementation Parity

Updated: 2026-06-22

## Purpose

This note records the implementation-parity checks for the full 4M neighborhood differential-abundance analyses. The manuscript should distinguish three related but non-identical analyses:

1. The project Python Milo-style workflow, which constructs stratified seed neighborhoods in the validated full 4M scVI latent space and fits donor-level logit-fraction age models.
2. An edgeR quasi-likelihood count-model parity workflow, which tests the same exported Python neighborhoods as neighborhood-by-donor counts.
3. An official MiloR subset sensitivity workflow, which constructs independent MiloR neighborhoods on a stratified subset of the same full 4M scVI latent space.

The goal is not to make all three workflows numerically identical. The goal is to test whether the biological direction survives changes in the statistical implementation and neighborhood construction.

## Current Status

- Python full 4M lineage DA is complete: 5,613 / 20,000 neighborhoods pass age FDR < 0.10.
- Python matched FLEX v2/device lineage DA is complete: 1 / 20,000 neighborhood passes age FDR < 0.10, labeled Early iOSN.
- Age-bin robustness is complete and directionally supports the matched Early iOSN result.
- edgeR count-model parity is complete for all-donor and matched lineage neighborhoods.
- Official MiloR subset sensitivity is complete using a user-level micromamba environment on `mia`.

## edgeR Count-Model Parity Results

The edgeR parity workflow exports the same Python neighborhoods as a neighborhood-by-donor count matrix and tests age with an edgeR quasi-likelihood model. This preserves the Python neighborhood definitions while changing the statistical core from donor logit-fraction OLS to negative-binomial count modeling.

| Run | Donors | Python FDR < 0.10 | edgeR FDR < 0.10 | Significant overlap | Signed-effect Spearman | Top-100 direction agreement |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Full 4M lineage | 187 | 5,613 | 4,758 | 1,375 | 0.916 | 1.00 among Python top 100; 0.92 among edgeR top 100 |
| Matched 4M lineage | 27 | 1 | 468 | 1 | 0.924 | 1.00 among Python top 100; 1.00 among edgeR top 100 |

The matched parity result is the strongest direct support for the current narrow claim: the single Python-significant matched Early iOSN neighborhood is also edgeR-significant, and signed age effects are highly concordant across all matched neighborhoods. The edgeR model calls more matched neighborhoods significant, so the manuscript should not use edgeR as a stricter filter. Instead, it should use edgeR as evidence that the direction of the Python donor-level result is not an artifact of the OLS/logit-fraction statistical core.

## Official MiloR Subset Results

Official MiloR was installed on `mia` in a self-contained micromamba environment because source installation against system R 4.6 failed on missing `fontconfig`/`freetype` development headers. The micromamba environment provides `bioconductor-milor` and its system-library dependencies without root access.

Official MiloR was run on stratified subsets exported from the same validated full 4M reduced scVI latent space:

| Run | Cells | Donors | MiloR neighborhoods | Spatial FDR < 0.10 | Covariates |
| --- | ---: | ---: | ---: | ---: | --- |
| Full 4M lineage subset | 100,000 | 187 | 4,579 | 3,447 | sex, chemistry, collection method, log total cells |
| Matched 4M lineage subset | 75,000 | 27 | 3,419 | 627 | sex, log total cells |

The official MiloR subsets confirm that the full 4M scVI lineage space contains strong age-associated differential-abundance structure. However, because MiloR constructs independent neighborhoods, its neighborhood rows are not one-to-one comparable with the Python seed neighborhoods. The dominant official MiloR subset signals are positive age-associated HBC/suprabasal/sustentacular neighborhoods. Negative significant MiloR neighborhoods are present, but in the matched subset they are mostly quiescent HBC and olfactory sustentacular neighborhoods; only one negative matched neighborhood is labeled Early mature mOSN, and Early iOSN is not a dominant matched official-MiloR signal.

This narrows the manuscript interpretation. Official MiloR supports the presence of age-associated latent-neighborhood structure but does not independently reproduce the exact matched Early iOSN neighborhood as the primary official-MiloR finding. The Early iOSN claim should therefore remain a narrow Python-neighborhood result supported by exact-neighborhood edgeR parity, age-bin directionality, and curated program enrichment, not a broad official-MiloR discovery.

## Rationale

MiloR performs graph-neighborhood differential abundance testing with negative-binomial generalized linear models or mixed models. The current project workflow intentionally differs in two ways: it uses preselected stratified seed neighborhoods to make full-scale runs reproducible and memory-bounded, and it fits donor-level logit-fraction models with donor covariates. This makes the Python workflow transparent and scalable for the full 4,028,275-cell atlas, but it must be presented as "Milo-style" rather than as official MiloR.

The edgeR parity workflow is therefore useful because it keeps the same Python neighborhoods but swaps the statistical core to a negative-binomial quasi-likelihood count model. The official MiloR subset workflow is useful because it swaps both the statistical core and the neighborhood constructor, while remaining computationally tractable.

## Acceptance Criteria

- edgeR parity produces neighborhood-by-donor count-model results for all-donor lineage and matched lineage runs. Complete.
- official MiloR subset sensitivity runs on at least one lineage subset from the full 4M scVI latent space. Complete for both full and matched lineage subsets.
- The final results note reports whether direction, top cell-state themes, and matched Early iOSN support are concordant, narrowed, or contradicted. Complete: official MiloR narrows rather than fully reproduces the matched Early iOSN interpretation.
- Manuscript text uses "Milo-style" for the Python workflow unless official MiloR results become the main analysis. Complete: official MiloR remains a sensitivity, not the main analysis.

## Interpretation Rules

- If edgeR parity and official MiloR subset both agree directionally with the Python lineage result, the neighborhood analysis can be promoted as a supported secondary mechanistic layer.
- If edgeR agrees but MiloR subset is weaker, keep the main-text claim conservative and place MiloR in Extended Data or Supplement. This is the current outcome.
- If either parity layer contradicts the matched Early iOSN interpretation, do not promote the neighborhood analysis beyond exploratory status.

## References

- Dann et al., "Differential abundance testing on single-cell data using k-nearest neighbor graphs," Nature Biotechnology, 2022. https://doi.org/10.1038/s41587-021-01033-z
- Bioconductor `miloR` package page. https://bioconductor.org/packages/miloR/
