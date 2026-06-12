# ORA Current Tasks and Milestone Tracker

Updated: 2026-06-12

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
- [x] Push/merge completed milestones to the private GitHub repo.

## Phase 1: External Validation

- [x] Add structured external dataset registry fields for independent olfactory aging or NDD datasets.
- [x] Add adapter CLI to summarize configured external datasets and required files.
- [x] Add external feature-matrix import contract for donor-level validation files.
- [x] Add published gene-list registry for aging, NDD, olfactory regeneration, and inflammation signatures.
- [x] Add gene-list validation command that reports coverage against Gateway genes.
- [ ] Add module scoring support for imported published gene lists.
- [ ] Test whether ORA-associated composition/module features replicate outside Gateway.
- [ ] Compare Gateway ORA signatures against external aging or disease contrasts.
- [ ] Report external validation status in the MVP report.

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

## Phase 3: Sensitivity Analyses

- [x] Re-run ORA by FLEX chemistry version.
- [x] Re-run ORA by collection method, especially brush vs Gateway device.
- [ ] Re-run targeted and adjusted DE by chemistry and collection method strata where sample size allows.
- [x] Add minimum cell-count threshold sweeps for donor/cell-state features.
- [x] Exclude low-cell donors and low-cell cell states in sensitivity runs.
- [x] Add healthy-only subcohort checks.
- [ ] Add bootstrap or leave-site/leave-chemistry-out sensitivity summaries.
- [ ] Add sensitivity tables and plots to the report.

## Phase 4: Improve ORA Modeling

- [ ] Add age calibration for ORA predictions.
- [ ] Add repeated donor-level CV.
- [ ] Add simpler interpretable models, such as ridge/lasso/linear baseline.
- [ ] Add feature stability selection across repeated CV.
- [ ] Add confidence intervals for MAE, Spearman r, and ORAA summaries.
- [ ] Compare composition-only vs module-augmented models formally.
- [ ] Add permutation/null-model tests.
- [ ] Add calibrated ORA plots and confidence intervals to the report.

## Phase 5: Interpret NDD Projection Carefully

- [ ] Add explicit warning/report note that AD/PD sample size is 5 donors each.
- [ ] Add sensitivity checks for negative ORAA.
- [ ] Separate disease biology from device, chemistry, and sample-composition effects.
- [ ] Compare NDD projection using composition-only vs module-augmented features.
- [ ] Add donor-level NDD projection appendix table.
- [ ] Add projection uncertainty or bootstrap intervals.
- [ ] Validate negative ORAA pattern in independent NDD olfactory datasets if available.

## Phase 6: Trajectory and Neighborhood Work

- [ ] Confirm whether original scANVI/scVI latent embeddings can be obtained outside CELLxGENE export.
- [ ] If unavailable, add workflow to recompute latent representations.
- [ ] Add pseudotime/lineage-density workflow using an appropriate latent space.
- [ ] Add Milo neighborhood analysis external-compute workflow.
- [ ] Add cNMF program discovery external-compute workflow.
- [ ] Add lineage bottleneck and density summaries.
- [ ] Add trajectory/neighborhood report sections only after latent-space validation.

## Phase 7: Manuscript/Package Hardening

- [ ] Add a reproducible `make all` or Snakemake profile for full real-data reruns.
- [ ] Add runtime/memory notes for real Gateway commands.
- [ ] Add data provenance table for every generated output.
- [ ] Add release checklist for figures/tables.
- [ ] Add CI-friendly tests that avoid requiring real Gateway data.
- [ ] Add method text drafts for ORA, module scoring, pseudobulk, and projection.

## Current Priority

1. External validation scaffold and published gene-list registry.
2. Genome-wide pseudobulk export hooks.
3. Sensitivity analyses for chemistry/collection method.
4. ORA model calibration and repeated CV.

## Progress Log

- 2026-06-12: Added this living task tracker.
- 2026-06-12: Added external validation registry, donor-level feature contract, published gene-list coverage command, and tests. Initial curated gene lists resolve 36/36 genes in Gateway; external datasets remain file-pending.
- 2026-06-12: Added genome-wide pseudobulk export CLI with gene chunking, sparse CSR aggregation, R DE hook templates, minimum group filters, and toy AnnData smoke coverage. Full Gateway export completed: 4,028,275 cells, 18,127 genes, 15,193 total groups, and 6,509 DE-ready groups. Local R-side DE remains pending because `Rscript` is not installed in this environment.
- 2026-06-12: Added genome-wide pseudobulk QC summaries. The real export has matching matrix/metadata IDs, 12,358 median detected genes per pseudobulk group, 0.8037 median gene detection fraction across groups, and 0.9994 matrix-to-metadata total-count ratio.
- 2026-06-12: Added ORA sensitivity reruns for chemistry, collection method, healthy-only, and minimum-cell thresholds. Ten scenarios ran successfully; random forest MAE was best in FLEX v2 donors (8.89, n=27) and weakest in brush-only donors (17.34, n=94), indicating collection method/chemistry sensitivity needs careful interpretation.
