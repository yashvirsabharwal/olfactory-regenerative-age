# ORA Current Tasks and Milestone Tracker

Updated: 2026-06-14

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
- [x] Add module scoring support for imported published gene lists.
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
- [ ] Add limma-voom and DESeq2 parity runs or documented adapter examples.

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
- [x] Add repeated donor-level CV.
- [ ] Add simpler interpretable models, such as ridge/lasso/linear baseline.
- [x] Add feature stability selection across repeated CV.
- [x] Add confidence intervals for MAE, Spearman r, and ORAA summaries.
- [ ] Compare composition-only vs module-augmented models formally.
- [ ] Add permutation/null-model tests.
- [ ] Add calibrated ORA plots and confidence intervals to the report.

## Phase 5: Interpret NDD Projection Carefully

- [x] Add explicit warning/report note that AD/PD sample size is 5 donors each.
- [ ] Add sensitivity checks for negative ORAA.
- [x] Separate disease biology from device, chemistry, and sample-composition effects.
- [ ] Compare NDD projection using composition-only vs module-augmented features.
- [ ] Add donor-level NDD projection appendix table.
- [x] Add projection uncertainty or bootstrap intervals.
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
