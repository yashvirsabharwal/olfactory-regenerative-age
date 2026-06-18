# ORA Methods Draft

Updated: 2026-06-18

## Cohort Construction

Gateway H5AD metadata are inspected in backed mode. Donor and sample columns are resolved from `configs/gateway.yaml`, disease labels are normalized into healthy, AD, PD, and other/unknown groups, and ORA training is restricted to healthy donors with valid age. AD/PD donors are held out from all age-model training and used only for frozen-model projection.

## Donor-Level Features

Fine-cell-state counts are aggregated by donor and transformed into composition features using proportions, centered log-ratio features, and curated lineage ratios. Technical variables such as chemistry, collection method, site, and yield are retained for diagnostics and sensitivity analyses, not used as primary biological ORA features.

## Module Scoring

Curated gene modules are scored by chunked average log1p expression over available genes, then summarized at donor and cell-state levels. Gene coverage is reported for every module before interpretation.

## External Validation Checks

For GSE184117, GEO sample metadata are parsed into donor age, olfaction status, sample class, and raw 10x file prefixes. Public raw matrices are used for sample-level module scoring, marker-only coarse composition scoring, marker-reference mapped AnnData/donor-feature generation, and Gateway scANVI/scArches query mapping. The scANVI query mapping intersects genes with the scaled Gateway reference, predicts Gateway fine cell-state labels, records label confidence and entropy, and exports ORA-compatible `prop__` and `clr__` donor features. Mapped features are compared with Gateway age-associated feature directions. These analyses remain claim-gated because the usable biopsy contrast is n=3 normosmic versus n=3 presbyosmic samples.

## ORA Modeling

ORA models predict chronological age from donor-level biological features in healthy donors only. Cross-validation is donor-level. Repeated CV, shuffled-age null models, nested tuning, calibration diagnostics, and residual diagnostics are used to distinguish robust signal from model instability.

## NDD Projection

Frozen healthy-trained ORA models are projected onto AD/PD donors. ORAA is interpreted only as a relative projected tissue-state deviation, with explicit caveats for n=5 AD, n=5 PD, and FLEX v2/device confounding.

## Pseudobulk DE

Targeted and genome-wide pseudobulk analyses aggregate counts by donor/sample/cell state. Genome-wide DE uses edgeR quasi-likelihood models with age, sex, chemistry, collection method, and site where estimable. limma-voom is run as a cross-method parity analysis on both the all-donor and matched FLEX v2/device subsets. DESeq2 is documented as deferred for the genome-wide fine-cell-state grid because donor balance and confounding, rather than absence of a third count model, are the limiting interpretability issues. Results require donor-balance, sentinel-gene, matched-subset, and edgeR/limma parity audits before biological interpretation.

## Latent-Space Validation

The Gateway CELLxGENE export is audited for `.obsm` embeddings before trajectory or neighborhood analysis. The public H5AD currently exposes `X_umap` but not a non-UMAP latent representation. To recover a latent atlas, the full Gateway H5AD is reduced in contiguous row chunks to a 3,003-gene HVG/marker feature set, concatenated on disk, and used to train scVI across all 4,028,275 cells with `sample_id` as the batch key and `flex_version`, `device_guided`, and `sex` as categorical covariates. Stratified 250,000-cell seed runs and a lineage-focused 100,000-cell run provide sensitivity anchors. Validation checks embedding dimensions, finite values, metadata coverage, nearest-neighbor label purity, technical mixing, and marker continuity. Full-scale Milo-style neighborhood analyses use the all-cell reduced scVI latent atlas with 20,000 stratified seed neighborhoods, 100 nearest cells per neighborhood, donor-level logit-fraction models, and adjustment for age, sex, chemistry, collection method, and donor yield. These neighborhood results remain mechanistic and claim-gated until matched technical sensitivity, marker/program annotation, and implementation-parity checks mature; pseudotime and cNMF remain deferred.
