# ORA Current Tasks and Milestone Tracker

Updated: 2026-06-18

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
- [x] Push completed milestone branch to the private GitHub repo.
- [x] Fast-forward or merge completed milestones into the private repo default branch.

## Phase 1: External Validation

- [x] Add structured external dataset registry fields for independent olfactory aging or NDD datasets.
- [x] Add adapter CLI to summarize configured external datasets and required files.
- [x] Add external feature-matrix import contract for donor-level validation files.
- [x] Add published gene-list registry for aging, NDD, olfactory regeneration, and inflammation signatures.
- [x] Add gene-list validation command that reports coverage against Gateway genes.
- [x] Add module scoring support for imported published gene lists.
- [ ] Add real external dataset download/import notes for each registry entry.
- [ ] Add raw single-cell adapter for external AnnData/H5AD inputs.
- [x] Add donor-feature adapter for external CSV/TSV feature matrices.
- [x] Add feature harmonization report for missing/renamed cell states and modules.
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
- [x] Add matched healthy reference DE for FLEX v2/device-only AD/PD comparisons.
- [x] Add per-cell-state donor balance diagnostics before interpreting DE.
- [x] Add a sex-linked and mitochondrial/ribosomal hit audit table.
- [x] Add limma-voom parity runs for all-donor and matched FLEX v2/device genome-wide DE.
- [x] Add DESeq2 parity example or documented reason if the design is too underpowered for stable DESeq2 fits.

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
- [x] Add sensitivity tables and plots to the report.

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
- [x] Add final model-card table describing features, exclusions, cohorts, and limitations.
- [x] Add calibrated ORA plots and confidence intervals to the report.

## Phase 5: Interpret NDD Projection Carefully

- [x] Add explicit warning/report note that AD/PD sample size is 5 donors each.
- [x] Add sensitivity checks for negative ORAA.
- [x] Separate disease biology from device, chemistry, and sample-composition effects.
- [x] Compare NDD projection using composition-only vs module-augmented features.
- [x] Add donor-level NDD projection appendix table.
- [x] Add disease projection residual diagnostics by age, sex, chemistry, collection method, site, and yield.
- [x] Add matched healthy reference plots for FLEX v2/device-only donors.
- [x] Add disease projection permutation test with labels shuffled within compatible strata where possible.
- [x] Add projection uncertainty or bootstrap intervals.
- [ ] Validate negative ORAA pattern in independent NDD olfactory datasets if available.

## Phase 6: Trajectory and Neighborhood Work

- [x] Confirm whether original scANVI/scVI latent embeddings can be obtained from the local H5AD or CELLxGENE portal metadata.
- [x] If unavailable, add workflow to recompute latent representations.
- [x] Add preprocessing plan for recomputed latent space, including HVG selection, covariates, and batch keys.
- [ ] Add pseudotime/lineage-density workflow using an appropriate latent space.
- [x] Add Milo neighborhood analysis external-compute workflow.
- [ ] Add cNMF program discovery external-compute workflow.
- [ ] Add lineage bottleneck and density summaries.
- [x] Validate that neighborhood findings are not UMAP-only artifacts.
- [x] Add full 4M Milo-style neighborhood report section after latent-space validation.

## Phase 7: Manuscript/Package Hardening

- [ ] Add a reproducible `make all` or Snakemake profile for full real-data reruns.
- [ ] Add runtime/memory notes for real Gateway commands.
- [x] Add data provenance table for every generated output.
- [ ] Add release checklist for figures/tables.
- [ ] Add CI-friendly tests that avoid requiring real Gateway data.
- [x] Add method text drafts for ORA, module scoring, pseudobulk, and projection.
- [x] Add paper claim ledger separating supported findings, exploratory observations, and deferred claims.
- [x] Add figure plan for main figures, extended data, and supplement tables.
- [x] Add limitations section draft focused on NDD sample size, chemistry/device confounding, and external validation.
- [x] Add reproducibility appendix with exact commands, software versions, and generated-output checksums.

## Phase 8: Paper Story and Novelty

- [x] Define primary claim: olfactory epithelial composition encodes a regenerative aging axis in healthy donors.
- [x] Define secondary claim: NDD projections are exploratory and hypothesis-generating only.
- [x] Prioritize interpretable cell-state signals over black-box prediction performance.
- [x] Map top ORA features to regenerative lineage biology, immune aging, glandular states, and neuronal maturation.
- [x] Compare ORA conceptually against generic tissue aging clocks and olfactory dysfunction biomarker literature.
- [x] Draft title, abstract, and one-page significance statement.
- [ ] Decide target venue tier after sensitivity/external-validation results mature.

## Needs From User / External Blockers

- [ ] Restore GitHub authentication so local commits can be pushed to the private repository.
- [ ] Provide or approve downloads for independent olfactory aging/NDD datasets when available.
- [ ] Provide any preferred target journal, paper style, or claim aggressiveness constraints.
- [ ] Provide paper PDF/supplement updates if a newer Gateway version appears.

## Current Priority

1. Resolve `GSE184117` cell labels or add defensible reference mapping, then convert raw 10x matrices into Gateway-compatible donor composition features.
2. Recover or recompute a non-UMAP latent space for trajectory/neighborhood analyses.
3. Resolve independent AD/PD olfactory validation data or keep NDD strictly exploratory.
4. Keep DESeq2 deferred unless independent disease donor counts improve or reviewers request a targeted subset.

## Detailed Remaining Work Queue

### A. External Validation Tasks

- [ ] For each external dataset candidate, add source URL, accession, download command, license/access notes, expected files, and donor/sample count.
- [x] Add first concrete external validation source notes for `GSE184117` and `GSE151973`.
- [x] Parse `GSE184117` GEO series-matrix metadata into sample/donor age, olfaction status, sample class, and raw 10x file prefixes.
- [x] Add sample-level raw 10x module scoring for `GSE184117` as a descriptive external sanity check.
- [x] Add marker-only coarse composition scoring for `GSE184117` as a descriptive external sanity check while public cell labels are unavailable.
- [x] Add marker-only concordance check comparing GSE184117 presbyosmia shifts with Gateway age-association directions.
- [x] Build raw external AnnData/H5AD adapter that resolves donor ID, age, disease, sex, chemistry/batch, cell labels, counts, and gene symbols.
- [x] Build external donor-feature matrix adapter for studies that only provide summarized composition/module tables.
- [x] Add external feature harmonization report covering missing cell-state aliases, renamed labels, unavailable modules, and feature drop rates.
- [x] Add claim-gated external validation evidence ledger separating direct validation candidates, raw-adapter candidates, marker-only context, and blocked datasets.
- [x] Add external composition/module feature generation using the same biological-feature contract as Gateway.
- [x] Run mapped-feature direction concordance tests for age-associated Gateway features in external presbyosmia samples.
- [ ] Run ORA signature replication tests for age-associated features in larger external healthy donors.
- [ ] Run external NDD projection or contrast tests if independent AD/PD olfactory data are available.
- [x] Add external validation report section separating ready, blocked, failed, marker-only, and raw-adapter-ready datasets.

### B. Genome-Wide DE Sensitivity Tasks

- [x] Add donor-balance diagnostics per contrast and fine cell state before interpreting genome-wide DE.
- [x] Add sex-linked, mitochondrial, ribosomal, hemoglobin, and immunoglobulin audit flags to genome-wide DE summaries.
- [x] Re-run AD/PD genome-wide DE against matched healthy FLEX v2/device donors where donor counts are sufficient.
- [ ] Run sex-stratified or sex-balanced genome-wide DE sensitivity where sample size allows.
- [ ] Run chemistry-stratified genome-wide DE sensitivity where sample size allows.
- [ ] Run collection-method-stratified genome-wide DE sensitivity where sample size allows.
- [x] Add limma-voom parity run on the exported pseudobulk matrix.
- [x] Add matched FLEX v2/device limma-voom sensitivity run.
- [x] Add DESeq2 parity example or documented reason if the design is too underpowered for stable DESeq2 fits.
- [x] Add report table that separates robust non-sex-linked signals from likely confounded sentinel hits.

### C. ORA/NDD Sensitivity Tasks

- [x] Add NDD matched FLEX v2/device healthy reference plot.
- [ ] Add matched healthy FLEX v2/device-only ORA training/projection sensitivity.
- [ ] Add leave-collection-method-out ORA checks.
- [ ] Add leave-chemistry-out or leave-site-out ORA checks where strata contain enough donors.
- [x] Add donor-yield sensitivity excluding top and bottom cell-count donors.
- [x] Add NDD disease-label permutation test within compatible FLEX v2/device strata if exchangeability is defensible.
- [ ] Add NDD projection figure captions and caveats directly into manuscript-ready notes.

### D. Modeling/Interpretability Tasks

- [x] Add final model-card table with model family, feature set, training cohort, exclusions, covariates, MAE, Spearman r, calibration slope, and limitations.
- [ ] Add calibrated ORA plots and confidence interval ribbons to the report.
- [x] Add histogram gradient boosting benchmark as a leakage-safe sklearn tabular model family.
- [ ] Decide whether TabPFN/AutoGluon is worth an exploratory run; if run, label it as non-primary because donor n is only 187.
- [x] Add top-feature biological interpretation table mapping ORA features to progenitor, neuronal maturation, immune, glandular, and injury/regeneration themes.

### E. Trajectory/Neighborhood Tasks

- [x] Search local Gateway H5AD and CELLxGENE portal metadata for original scVI/scANVI latent embeddings.
- [x] If latent embeddings are unavailable, add reproducible scVI recomputation workflow with HVG, batch covariates, chemistry, and donor/sample handling.
- [x] Install optional latent dependencies and run a pilot scVI recomputation on a bounded subset.
- [x] Validate pilot latent space with marker continuity and donor/chemistry/collection-method mixing diagnostics.
- [x] Add scaled 250k-cell and lineage-focused scVI Make targets with stratified sampling.
- [x] Run scaled 250k-cell stratified scVI and validate marker continuity, label purity, and metadata mixing.
- [x] Run lineage-focused basal/progenitor/neural scVI and validate before pseudotime.
- [x] Train a full 4,028,275-cell reduced Gateway scVI model on 3,003 HVG/marker genes.
- [x] Run full 4M Milo-style neighborhood DA for broad all-cell, lineage-focused, and secretory-focused neighborhoods.
- [x] Run matched FLEX v2/device full 4M Milo-style sensitivity.
- [x] Annotate significant neighborhoods by ORA feature themes and claim gates.
- [ ] Add gene-level marker/program enrichment for significant neighborhoods.
- [ ] Add official MiloR parity or a documented implementation-parity rationale.
- [ ] Add pseudotime or lineage-density workflow only after latent space is validated.
- [x] Add first-pass Milo-style neighborhood workflow after latent space validation.
- [x] Run publication-scale Milo-style workflow on the full 4M reduced scVI atlas.
- [x] Run matched/sensitivity Milo-style neighborhood workflow before biological claims.
- [ ] Add cNMF program workflow only after latent space is validated.
- [ ] Add negative-control check showing conclusions are not UMAP-only artifacts.

### F. Manuscript/Package Hardening Tasks

- [ ] Add `Makefile` or Snakemake profile for full real-data reruns.
- [x] Add command manifest with every real-data command, input, output, runtime, memory note, and software environment.
- [x] Add output provenance/checksum table for generated report tables and figures.
- [x] Add CI-friendly smoke tests that avoid real Gateway data.
- [x] Draft methods sections for cohort construction, ORA modeling, module scoring, pseudobulk DE, and NDD projection.
- [x] Draft limitations section covering external validation, NDD n=5, chemistry/device confounding, sex imbalance, and missing latent embeddings.
- [x] Draft claim ledger separating supported findings, exploratory findings, negative benchmarks, and deferred claims.
- [x] Draft figure plan for main figures, extended data, and supplement tables.

### G. Paper Story/Novelty Tasks

- [x] Define primary claim in one sentence: olfactory epithelial composition encodes a regenerative aging axis in healthy donors.
- [x] Define secondary claim in one sentence: NDD projections are exploratory and hypothesis-generating.
- [x] Decide whether the manuscript is best framed as resource/methods, atlas reanalysis, aging biology, or NDD hypothesis paper.
- [x] Compare ORA conceptually against tissue clocks, olfactory dysfunction aging literature, and regenerative epithelial-state papers.
- [x] Draft title, abstract, and one-page significance statement after external validation status is clear.

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
- 2026-06-15: Added NDD projection residual diagnostics by sex, age bin, chemistry, collection method, site, and cell-yield quartile, with single-donor strata explicitly flagged. The real augmented projection generated 288 diagnostic rows. Chemistry and collection method are fully confounded for disease donors because all AD/PD donors are FLEX v2/device; site is missing in the export. Among ok strata for the current display models, PD female and high-yield strata remain strongly negative, and AD male/70-79 strata remain negative, but many disease strata are too small for inference.
- 2026-06-15: Added a matched FLEX v2/device healthy-reference NDD ORAA figure to the MVP report and expanded the task tracker into a granular remaining-work queue. External validation research identified GEO `GSE184117` as the first actionable human olfactory aging/presbyosmia single-cell validation target and `GSE151973` as a bulk olfactory neuroepithelium marker-sanity reference; source notes and download commands are now tracked in `docs/external_validation_sources.md`.
- 2026-06-15: Added reproducibility hardening and manuscript-readiness summaries. The Makefile now prefers `.venv/bin/python`, ruff hygiene is clean, external feature-matrix validation and raw-archive inventory commands are available, `GSE184117_RAW.tar` was downloaded and inventoried as 21 10x-style files across 3 controls, 3 presbyosmic samples, and 1 culture sample, genome-wide DE audit/donor-balance/matched-feasibility tables are generated, NDD frozen-score label permutation is reported, an ORA model-card table is generated, output provenance now tracks 58 artifacts with no missing outputs, and draft methods/claim/figure/limitations/reproducibility documents are in `docs/`.
- 2026-06-15: Parsed `GSE184117` GEO series-matrix metadata and added sample-level raw 10x module scoring. The external validation registry now marks `oliva_2022` as raw-adapter-ready, with 6 usable biopsy samples (3 normosmic, 3 presbyosmic/hyposmic/anosmic) plus one excluded culture sample. Fourteen curated/published modules were scored directly from the raw 10x matrices; contrasts remain descriptive because n=3 versus n=3 and public cell labels are unavailable.
- 2026-06-15: Added histogram gradient boosting as a leakage-safe sklearn tabular model family and ran a targeted 10-repeat module-augmented candidate benchmark against XGBoost, CatBoost, and the boosted ensemble. CatBoost remained best by mean MAE (14.08), histogram gradient boosting was essentially tied by MAE (14.09) and had the highest candidate Spearman r (0.371), so it strengthens the benchmark set but does not change the main conclusion that accuracy gains are modest.
- 2026-06-15: Added calibrated ORA report figures and ran matched FLEX v2/device genome-wide edgeR sensitivity. The matched subset used 37 donors and 1,942 pseudobulk groups, tested 852,517 rows across 82 successful cell-state contrast models, and reduced total FDR-significant rows from 1,658 to 641. AD matched DE has 351 significant rows and 0 sex-linked sentinel rows; PD matched DE has 290 significant rows but still 168 sex-linked sentinel rows, reinforcing that PD genome-wide biology remains confounding-sensitive.
- 2026-06-15: Added genome-wide limma-voom parity for both all-donor and matched FLEX v2/device pseudobulk DE. All-donor limma-voom tested 862,707 rows and found 174 FDR-significant rows, versus 1,658 under edgeR; matched limma-voom tested 852,517 rows and found 127 significant rows, with 0 sex-linked sentinel rows in both AD and PD matched contrasts. This strengthens the DE section as a conservative cross-method sensitivity analysis and keeps biological claims focused on robust, audited signals.
- 2026-06-15: Added marker-only GSE184117 coarse composition scoring over 11 olfactory/epithelial/immune marker panels plus an unassigned bin. The real archive produced 77 marker-coverage rows, 84 sample/panel composition rows, and 12 tiny-n marker contrasts. Presbyosmia samples were descriptively higher for goblet/secretory, mature OSN, cycling, activated HBC, and quiescent HBC marker fractions, but all contrasts remain `marker_only_small_n` because n=3 versus n=3 and no public cell labels are available.
- 2026-06-15: Documented DESeq2 as a deferred genome-wide parity engine in `docs/deseq2_parity_decision.md`. The rationale is that edgeR QL plus limma-voom already provide cross-method sensitivity, while the main disease-DE limitation is donor balance/confounding rather than absence of a third count model. DESeq2 should be revisited for larger independent disease cohorts, reviewer-requested shortlisted hits, or a small predeclared cell-state/gene subset.
- 2026-06-15: Added a latent-space readiness audit. The local Gateway H5AD has only `X_umap` (4,028,275 cells x 2 dimensions), and the CELLxGENE collection API reports only one H5AD asset with portal embeddings `X_umap`. The status is therefore `latent_recompute_required`; trajectory, Milo, and cNMF remain blocked until a non-UMAP latent representation is recovered from authors or recomputed and validated.
- 2026-06-15: Added a guarded scVI latent recomputation workflow scaffold. `make latent-space-recompute-plan` writes dependency/scale feasibility and `docs/latent_recompute_workflow.md`; `scripts/run_scvi_latent.py` is the explicit training entry point after installing the `latent` optional dependency extra.
- 2026-06-15: Added a claim-gated external validation evidence ledger. The real run generated 6 rows: 4 configured source rows plus GSE184117 sample-module and marker-only composition sanity checks. GSE184117 is explicitly marked as raw-adapter candidate plus sanity-only generated evidence; GSE151973 remains marker-context-only bulk support.
- 2026-06-15: Added a top-feature biological interpretation table. The real augmented repeated-CV stability run now feeds `results/tables/ora_feature_interpretation.tsv`, which maps 30 stable ORA features to supporting/secretory epithelium, immune/inflammatory, stress/senescence, neuronal-lineage, regenerative/progenitor, and disease-module themes with per-row cautions.
- 2026-06-15: Added `docs/manuscript_story.md`, a web-checked manuscript framing note with primary/secondary claims, novelty, external-source context, publication-readiness verdict, draft title, abstract, and one-page significance language. Current framing is atlas reanalysis/resource plus healthy aging biology, not disease-biomarker or trajectory paper.
- 2026-06-16: Added marker-only external concordance against Gateway age associations. GSE184117 marker-panel presbyosmia shifts now map to 32 Gateway marker-feature rows; 13 are directionally concordant and 19 discordant, preserving the sanity-only interpretation because the external cohort is still 3 versus 3 without cell labels.
- 2026-06-16: Installed latent dependencies and completed a bounded scVI pilot on the Gateway H5AD. The run sampled 25,000 cells, selected 2,000 HVGs, trained 20 epochs with `sample_id` as the scVI batch key and `flex_version/device_guided/sex` covariates, and wrote `data/processed/gateway_scvi_pilot_25k.h5ad` locally. `results/tables/scvi_pilot_validation.tsv` confirms `X_scvi` is present (25,000 x 10) and finite, while sparse donor/sample/fine-cell-type levels keep trajectory, Milo, and cNMF claims gated until a larger or lineage-focused validation run passes marker-continuity and mixing diagnostics.
- 2026-06-16: Added a full manuscript framework and polished main-text figure builder. `docs/manuscript_framework.md` now defines the exact section structure, figure mapping, claim hierarchy, and publication standing. `make manuscript-figures` writes six manuscript-oriented PNG/PDF figures covering cohort/design, aging composition, ORA modeling, feature biology, external/NDD guardrails, and DE/latent readiness.
- 2026-06-16: Drafted the full LaTeX manuscript in `manuscript/main.tex` with DOI-backed citations in `manuscript/references.bib`, six main figure callouts, methods, limitations, and data/code availability. Citation keys resolve and brace-level LaTeX sanity checks pass; PDF compilation is pending installation of `latexmk` or `pdflatex`/`bibtex` on the local machine.
- 2026-06-16: Upgraded the latent recomputation path for the max-novelty manuscript plan. `scripts/run_scvi_latent.py` now supports deterministic stratified sampling across biological/technical metadata, lineage-focused regex filters, and marker-preserving HVG selection through Gateway `feature_name` symbols. `make scvi-scaled-250k` and `make scvi-lineage-basal-neural` are registered with validation/provenance targets. The regenerated 25k pilot now retains 2,009 variables including canonical markers; validation reports strong fine/coarse label purity, acceptable FLEX/device/condition/sex mixing, and interpretable basal, OSN, and sustentacular marker enrichment while keeping progenitor/immune and all trajectory/Milo/cNMF claims gated for scaled/lineage validation.
- 2026-06-16: Completed the scaled latent studies. `make scvi-scaled-250k` trained a 250,000-cell stratified Gateway scVI atlas with 3,003 HVGs/markers in ~30 minutes and wrote `data/processed/gateway_scvi_stratified_250k.h5ad`; validation confirms finite 250,000 x 10 `X_scvi`, strong fine/coarse label purity, acceptable FLEX/device/condition/sex mixing, and interpretable basal, immature OSN, mature OSN, sustentacular, and immune marker continuity, with progenitor markers still limited. `make scvi-lineage-basal-neural` trained a 100,000-cell lineage model with 3,008 HVGs/markers in ~11 minutes and similarly passed label-purity and technical-mixing gates, but progenitor and immune panels remain limited. `make output-provenance` now reports 125 outputs with 0 missing.
- 2026-06-18: Added a first-pass Python Milo-style neighborhood DA workflow. `make milo` now tests donor-level age association across scVI neighborhoods from the 250k atlas, adjusting for sex, chemistry, and collection method. The first run tested 1,000 neighborhoods and found 0 age-associated neighborhoods at BH FDR < 0.10; the strongest nominal hit was quiescent HBC-enriched but not claim-ready.
- 2026-06-18: Added targeted `make milo-lineage` and `make milo-secretory` passes. Both tested 1,000 neighborhoods and found 0 age-associated neighborhoods at BH FDR < 0.10. The strongest nominal lineage neighborhoods were early mature mOSN, late iOSN, and cycling HBC with negative age coefficients; the strongest nominal secretory neighborhoods were serous, club, mucous gland, and proliferating secretory with negative age coefficients. These remain exploratory until matched sensitivity and/or full-model neighborhood transfer improve support.
- 2026-06-18: Promoted the next Milo-style step to full-data scale. `make milo-full-4m`, `make milo-full-4m-lineage`, and `make milo-full-4m-secretory` target `data/processed/gateway_scvi_full_4m_reduced.h5ad` with 20,000 stratified neighborhoods, 100 nearest cells per neighborhood, minimum 30 donors, and donor-yield adjustment. These are the intended publishable runs; smaller atlas runs remain method checks and sensitivity anchors.
- 2026-06-18: Completed the full 4M Milo-style run set on `mia`. The broad full run found 3,673/20,000 age-associated neighborhoods at FDR < 0.10, the lineage-focused run found 5,613/20,000, and the secretory-focused run found 285/20,000. The clearest biological signal is negative age association in immature olfactory neuronal/regenerative neighborhoods, with broad all-cell support for secretory/ciliated/immune neighborhood shifts. Added `docs/milo_full_4m_results.md` to capture interpretation and remaining claim gates.
- 2026-06-18: Completed matched FLEX v2/device full 4M Milo-style sensitivity. With 27 healthy matched donors, the all-cell run retained 10 FDR < 0.10 neighborhoods, mostly negative Naive CD8/T-cell neighborhoods plus one late iOSN neighborhood; the lineage-focused run retained one negative Early iOSN neighborhood; and the secretory-focused run retained 0 significant neighborhoods. This narrows the manuscript-ready mechanistic language to technically matched regenerative neuronal-lineage depletion and keeps secretory-only neighborhood DA exploratory.
- 2026-06-18: Added full 4M Milo-style theme annotation. `make milo-full-4m-annotation` writes `milo_full_4m_top_neighborhood_themes.tsv` and `milo_full_4m_theme_summary.tsv`, mapping significant neighborhoods to ORA biology themes and claim gates. The matched regenerative-neuronal theme is now explicitly separated from all-donor hypothesis-map themes and matched immune-support themes. Gene-level marker/program enrichment remains a separate open task because neighborhood membership was not emitted by the DA runner.
- 2026-06-18: Added optional neighborhood membership export to the Milo-style runner. Passing `--membership-out` now writes exact neighborhood-to-cell membership rows with original cell indices and obs names, making the next marker/program enrichment pass feasible on `mia` without redesigning the DA workflow.
- 2026-06-16: Added and ran the GSE184117 raw external AnnData adapter. `make external-gse184117-mapped` wrote `data/processed/gse184117_marker_mapped.h5ad` with 59,656 biopsy cells x 36,601 genes, `results/tables/external_10x_mapping_qc.tsv` for 6 usable samples, and `data/processed/gse184117_mapped_donor_features.tsv` with 25 ORA-compatible `prop__`, `clr__`, and `ratio__` features. Mapping QC reports 95.6-99.9% mapped cells and 43 marker genes present per sample. `make external-mapped-feature-concordance` generated 32 Gateway-mapped feature-direction rows: 16 concordant and 16 discordant, all explicitly `small_n_mapped`. The external evidence ledger now includes the mapped-feature candidate row, and output provenance reports 129 outputs with 0 missing.
