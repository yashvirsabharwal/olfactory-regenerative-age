# ORA Current Tasks and Milestone Tracker

Updated: 2026-06-15

This is the working task board for the olfactory-regenerative-age project. Checkboxes should be updated as implementation lands, tests pass, and real Gateway outputs are regenerated.

## Foundation Completed

- [x] Scaffold reproducible ORA repository.
- [x] Add Gateway paper-aligned defaults and metadata mappings.
- [x] Inspect Gateway H5AD in backed mode.
- [x] Build cohort manifest and cohort summary.
- [x] Aggregate donor/cell-state composition features.
- [x] Build composition-only ORA feature matrix.
- [x] Train healthy-only ORA age models with donor-level CV.
- [x] Generate MVP report with figures.
- [x] Score curated biology modules.
- [x] Build module-augmented ORA feature matrix and models.
- [x] Project frozen healthy-trained ORA models onto AD/PD donors.
- [x] Run targeted curated-gene pseudobulk DE.
- [x] Run covariate-adjusted targeted pseudobulk DE.
- [x] Keep generated data/results out of Git.
- [x] Create local milestone commits for completed work.
- [ ] Push/merge completed milestones to the private GitHub repo once GitHub auth is restored.

## Phase 1: External Validation

- [x] Add structured external dataset registry fields for independent olfactory aging or NDD datasets.
- [x] Add adapter CLI to summarize configured external datasets and required files.
- [x] Add external feature-matrix import contract for donor-level validation files.
- [x] Add published gene-list registry for aging, NDD, olfactory regeneration, and inflammation signatures.
- [x] Add gene-list validation command that reports coverage against Gateway genes.
- [x] Add module scoring support for imported published gene lists.
- [ ] Add real external dataset download/import notes for each registry entry.
- [ ] Add raw single-cell adapter for external AnnData/H5AD inputs.
- [ ] Add donor-feature adapter for external CSV/TSV feature matrices.
- [ ] Add feature harmonization report for missing/renamed cell states and modules.
- [ ] Test whether ORA-associated composition/module features replicate outside Gateway.
- [ ] Compare Gateway ORA signatures against external aging or disease contrasts.
- [x] Report external validation status in the MVP report.

## Phase 2: Genome-Wide Pseudobulk DE

- [x] Add genome-wide pseudobulk export that streams H5AD without dense matrix loading.
- [x] Support gene chunking and sparse CSR reads for all genes.
- [x] Emit wide count matrix and sample metadata formats compatible with edgeR/limma/DESeq2.
- [x] Add command templates for limma-voom, edgeR quasi-likelihood, and DESeq2.
- [x] Add design formulas for disease, age, sex, chemistry, collection method, and site where available.
- [x] Add minimum donor/cell count filters before export.
- [x] Add smoke tests on toy AnnData for genome-wide export.
- [x] Add report section summarizing genome-wide DE readiness and outputs.
- [x] Add genome-wide pseudobulk QC summaries for matrix alignment, detection rates, disease strata, and top variable genes.
- [x] Install local R/Bioconductor runtime with edgeR, limma, and DESeq2 via micromamba.
- [x] Run full genome-wide edgeR quasi-likelihood DE on the exported Gateway pseudobulk matrix.
- [x] Add compact genome-wide DE summary and top-hit report tables.
- [ ] Add sex-balanced or sex-stratified genome-wide DE sensitivity models where sample size allows.
- [ ] Add chemistry-stratified genome-wide DE sensitivity models where sample size allows.
- [ ] Add collection-method-stratified genome-wide DE sensitivity models where sample size allows.
- [ ] Add matched healthy reference DE for FLEX v2/device-only AD/PD comparisons.
- [ ] Add per-cell-state donor balance diagnostics before interpreting DE.
- [ ] Add a sex-linked and mitochondrial/ribosomal hit audit table.
- [ ] Add limma-voom and DESeq2 parity runs or documented adapter examples.

## Phase 3: Sensitivity Analyses

- [x] Re-run ORA by FLEX chemistry version.
- [x] Re-run ORA by collection method, especially brush vs Gateway device.
- [ ] Re-run targeted and adjusted DE by chemistry and collection method strata where sample size allows.
- [x] Add minimum cell-count threshold sweeps for donor/cell-state features.
- [x] Exclude low-cell donors and low-cell cell states in sensitivity runs.
- [x] Add healthy-only subcohort checks.
- [ ] Add bootstrap or leave-site/leave-chemistry-out sensitivity summaries.
- [ ] Add leave-collection-method-out model checks.
- [ ] Add matched healthy FLEX v2/device-only ORA training/projection checks.
- [ ] Add donor-yield sensitivity checks excluding top/bottom cell-count donors.
- [ ] Add sensitivity tables and plots to the report.

## Phase 4: Improve ORA Modeling

- [x] Add age calibration for ORA predictions.
- [x] Add repeated donor-level CV.
- [x] Add simpler interpretable models, such as ridge/lasso/linear baseline.
- [x] Benchmark stronger tree/tabular models beyond random forest.
- [x] Add optional native booster backends for XGBoost, LightGBM, and CatBoost with safe sklearn fallbacks.
- [x] Add feature stability selection across repeated CV.
- [x] Add confidence intervals for MAE, Spearman r, and ORAA summaries.
- [x] Compare composition-only vs module-augmented models formally.
- [x] Add permutation/null-model tests.
- [x] Add model residual diagnostics by sex, chemistry, collection method, site, and cell yield.
- [x] Add calibration slope/intercept and age-bin error summaries.
- [x] Add nested hyperparameter tuning for the top booster models without leaking donor folds.
- [x] Add leakage-safe out-of-fold stacking for top ORA model families.
- [ ] Evaluate pretrained/AutoML tabular options such as TabPFN or AutoGluon only as exploratory benchmarks, and only if they can be defended against overfitting at n=187 donors.
- [ ] Add final model-card table describing features, exclusions, cohorts, and limitations.
- [ ] Add calibrated ORA plots and confidence intervals to the report.

## Phase 5: Interpret NDD Projection Carefully

- [x] Add explicit warning/report note that AD/PD sample size is 5 donors each.
- [x] Add sensitivity checks for negative ORAA.
- [x] Separate disease biology from device, chemistry, and sample-composition effects.
- [x] Compare NDD projection using composition-only vs module-augmented features.
- [x] Add donor-level NDD projection appendix table.
- [ ] Add disease projection residual diagnostics by age, sex, chemistry, collection method, site, and yield.
- [ ] Add matched healthy reference plots for FLEX v2/device-only donors.
- [ ] Add disease projection permutation test with labels shuffled within compatible strata where possible.
- [x] Add projection uncertainty or bootstrap intervals.
- [ ] Validate negative ORAA pattern in independent NDD olfactory datasets if available.

## Phase 6: Trajectory and Neighborhood Work

- [ ] Confirm whether original scANVI/scVI latent embeddings can be obtained outside CELLxGENE export.
- [ ] If unavailable, add workflow to recompute latent representations.
- [ ] Add preprocessing plan for recomputed latent space, including HVG selection, covariates, and batch keys.
- [ ] Add pseudotime/lineage-density workflow using an appropriate latent space.
- [ ] Add Milo neighborhood analysis external-compute workflow.
- [ ] Add cNMF program discovery external-compute workflow.
- [ ] Add lineage bottleneck and density summaries.
- [ ] Validate that trajectory/neighborhood findings are not UMAP-only artifacts.
- [ ] Add trajectory/neighborhood report sections only after latent-space validation.

## Phase 7: Manuscript/Package Hardening

- [ ] Add a reproducible `make all` or Snakemake profile for full real-data reruns.
- [ ] Add runtime/memory notes for real Gateway commands.
- [ ] Add data provenance table for every generated output.
- [ ] Add release checklist for figures/tables.
- [ ] Add CI-friendly tests that avoid requiring real Gateway data.
- [ ] Add method text drafts for ORA, module scoring, pseudobulk, and projection.
- [ ] Add paper claim ledger separating supported findings, exploratory observations, and deferred claims.
- [ ] Add figure plan for main figures, extended data, and supplement tables.
- [ ] Add limitations section draft focused on NDD sample size, chemistry/device confounding, and external validation.
- [ ] Add reproducibility appendix with exact commands, software versions, and generated-output checksums.

## Phase 8: Paper Story and Novelty

- [ ] Define primary claim: olfactory epithelial composition encodes a regenerative aging axis in healthy donors.
- [ ] Define secondary claim: NDD projections are exploratory and hypothesis-generating only.
- [ ] Prioritize interpretable cell-state signals over black-box prediction performance.
- [ ] Map top ORA features to regenerative lineage biology, immune aging, glandular states, and neuronal maturation.
- [ ] Compare ORA conceptually against generic tissue aging clocks and olfactory dysfunction biomarker literature.
- [ ] Draft title, abstract, and one-page significance statement.
- [ ] Decide target venue tier after sensitivity/external-validation results mature.

## Needs From User / External Blockers

- [ ] Restore GitHub authentication so local commits can be pushed to the private repository.
- [ ] Provide or approve downloads for independent olfactory aging/NDD datasets when available.
- [ ] Provide any preferred target journal, paper style, or claim aggressiveness constraints.
- [ ] Provide paper PDF/supplement updates if a newer Gateway version appears.

## Current Priority

1. External validation with imported published signatures and any available independent olfactory/NDD datasets.
2. Sex/chemistry/collection-method sensitivity for genome-wide DE and NDD projection.
3. ORA model calibration, interpretable baselines, and formal composition-vs-module comparison.
4. Latent-space recovery or recomputation for trajectory/neighborhood analyses.

## Progress Log

- 2026-06-12: Added this living task tracker.
- 2026-06-12: Added external validation registry, donor-level feature contract, published gene-list coverage command, and tests. Initial curated gene lists resolve 36/36 genes in Gateway; external datasets remain file-pending.
- 2026-06-12: Added genome-wide pseudobulk export CLI with gene chunking, sparse CSR aggregation, R DE hook templates, minimum group filters, and toy AnnData smoke coverage. Full Gateway export completed: 4,028,275 cells, 18,127 genes, 15,193 total groups, and 6,509 DE-ready groups. Local R-side DE remains pending because `Rscript` is not installed in this environment.
- 2026-06-12: Added genome-wide pseudobulk QC summaries. The real export has matching matrix/metadata IDs, 12,358 median detected genes per pseudobulk group, 0.8037 median gene detection fraction across groups, and 0.9994 matrix-to-metadata total-count ratio.
- 2026-06-12: Added ORA sensitivity reruns for chemistry, collection method, healthy-only, and minimum-cell thresholds. Ten scenarios ran successfully; random forest MAE was best in FLEX v2 donors (8.89, n=27) and weakest in brush-only donors (17.34, n=94), indicating collection method/chemistry sensitivity needs careful interpretation.
- 2026-06-13: Added 10-repeat donor-level CV for ORA. Composition-only random forest mean MAE is 14.20 with empirical 95% interval 13.87-14.48, and mean Spearman r is 0.330 with interval 0.306-0.365.
- 2026-06-14: Added NDD projection bootstrap uncertainty and chemistry/device context. AD and PD each have 5 donors, all FLEX v2/device; matched healthy FLEX v2/device reference has 16 donors. Negative ORAA remains below matched healthy reference for random forest: AD mean -9.06, matched difference -8.54 (95% bootstrap interval -12.66 to -3.67); PD mean -12.17, matched difference -11.65 (-17.51 to -6.66).
- 2026-06-14: Installed a local micromamba R/Bioconductor environment (`.mamba/ora-r`) with edgeR 4.4.0, limma 3.62.1, DESeq2 1.46.0, data.table, optparse, and R.utils. Full genome-wide edgeR QL DE completed on the 18,127-gene export, producing 862,707 tested gene/cell-state/contrast rows across 83 successful models. At FDR < 0.05, AD vs healthy has 819 significant rows across 28 cell states and PD vs healthy has 839 significant rows across 38 cell states. PD top hits are dominated by Y-linked sentinel genes because the PD cohort is 1 male/4 female versus healthy 77 male/113 female/2 unknown; report now includes a non-sex-linked top-hit table for transparent triage.
- 2026-06-14: Connected published external gene lists to the module-scoring engine and MVP report. Full Gateway scoring completed for four published-list modules (aging, olfactory regeneration, neurodegeneration risk, neuroinflammation), generating 60,772 grouped score rows and 202-donor features. All 36 requested genes resolve in Gateway. External validation registry tracks three candidate datasets, but none are feature-ready until expression/metadata or donor feature matrices are supplied.
- 2026-06-14: Added ridge and lasso as interpretable ORA baselines across standard CV, repeated CV, feature importance, NDD projection, uncertainty summaries, and report figures. In 10-repeat composition-only CV, ridge has MAE 14.60 (95% interval 14.11-15.12) and Spearman r 0.302, lasso has MAE 15.12 and Spearman r 0.217, elastic net has MAE 14.75, and random forest remains best overall with MAE 14.20 and Spearman r 0.330.
- 2026-06-14: Added ORA calibration and residual diagnostics. Full composition-only score diagnostics now output calibration summaries, calibrated scores, age-bin errors, and residual strata by sex, chemistry, collection method, race/ethnicity, site, and cell-yield quartile. Random forest remains best by raw and recalibrated MAE (13.85 and 13.84), but all non-null models are under-dispersed with ORA-on-age calibration slopes around 0.11-0.18; young donors are overpredicted and old donors are underpredicted, so ORA should be interpreted as a relative aging/regenerative-state axis rather than an absolute age clock.
- 2026-06-14: Added Extra Trees, gradient boosting, and a conservative tree-ensemble average to test whether random forest was the practical ceiling. In the sklearn-only 10-repeat donor-level CV pass, random forest was best by MAE (14.20), while the tree ensemble was close by MAE (14.27) and slightly stronger on RMSE/R2/Spearman (RMSE 17.64, R2 0.107, Spearman 0.332). Gradient boosting was competitive but less stable by MAE (14.48); this result was later superseded by native booster benchmarking.
- 2026-06-14: Added optional XGBoost, LightGBM, CatBoost, and boosted-ensemble model families. XGBoost is now the top composition-only model in 20-repeat donor-level CV (mean MAE 14.15; 95% empirical interval 13.74-14.70; mean Spearman r 0.354), with random forest close behind (MAE 14.21). Module-augmented repeated CV is slightly better overall: CatBoost has the best mean MAE (14.08; 13.60-14.57), followed by augmented XGBoost (14.12). The module feature set improves MAE for most model families, but the gains are modest and the intervals overlap, so the current model claim should emphasize robust benchmarking rather than a decisive algorithmic breakthrough.
- 2026-06-14: Added formal composition-vs-module comparison outputs (`ora_feature_set_model_comparison.tsv` and `ora_feature_set_model_deltas.tsv`) and integrated the ranked comparison into the MVP report. NDD projection and bootstrap uncertainty were regenerated with the expanded model set; negative AD/PD ORAA remains present across model families, but the interpretation remains exploratory because AD and PD each have only 5 donors and share FLEX v2/device chemistry/collection context.
- 2026-06-15: Added donor-level shuffled-age permutation/null testing with model-subset support, a CLI, unit tests, config outputs, and MVP report integration. A 50-permutation, 2-repeat augmented-feature null run completed for random forest, XGBoost, CatBoost, and the boosted ensemble. All four observed models beat every shuffled-label permutation for MAE and Spearman r, giving empirical p=0.0196 at this permutation depth. CatBoost observed MAE is 14.08 versus shuffled null mean 15.99; XGBoost observed MAE is 14.12 versus null mean 16.56; boosted ensemble observed MAE is 14.16 versus null mean 16.37; random forest observed MAE is 14.19 versus null mean 16.02.
- 2026-06-15: Added leakage-safe nested hyperparameter tuning for top booster models with inner donor-level CV inside each outer fold, CLI outputs, unit tests, and MVP report integration. Native XGBoost/LightGBM required restoring the local `libomp.dylib` runtime link from `.mamba/ora-r`. A 5-repeat augmented-feature nested tuning run showed modest gains: tuned XGBoost mean MAE 14.04 (13.65-14.41) versus untuned repeated-CV mean 14.12, while tuned CatBoost mean MAE 14.08 (13.79-14.52) was essentially unchanged from untuned 14.08. Current conclusion: nested tuning helps a little for XGBoost but does not move the project anywhere near zero MAE.
- 2026-06-15: Added leakage-safe out-of-fold stacking across ridge, random forest, XGBoost, and CatBoost with inner donor-fold meta-training, CLI outputs, unit tests, and MVP report integration. A 5-repeat augmented-feature run produced stacked-ensemble mean MAE 14.46 (13.98-14.88) and mean Spearman r 0.313, which is worse than tuned XGBoost MAE 14.04 and tuned CatBoost MAE 14.08. The ensemble weights were unstable across folds, so stacking is recorded as a negative benchmark rather than the preferred ORA model.
- 2026-06-15: Added NDD projection feature-set sensitivity and donor appendix exports. The workflow reruns frozen ORA projection with composition-only and module-augmented features, emits a 24-row model/disease comparison, and writes a 240-row AD/PD donor appendix. Excluding the null model, all 11 model families show negative mean ORAA for both AD and PD under both feature sets; the largest feature-set shift was CatBoost AD, moving from -8.93 to -6.13 years. This supports directional stability of the exploratory NDD projection, but the sample size remains 5 AD and 5 PD donors.
