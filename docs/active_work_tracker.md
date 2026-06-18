# ORA Active Work Tracker

Updated: 2026-06-18

This is the live tracker for the post-preprint upgrade path. Mark tasks here as they move from planned to implemented, validated, and manuscript-ready.

## 0. Repository Hygiene

- [x] Scan tracked and ignored repository structure.
- [x] Remove generated Python/ruff cache directories.
- [x] Update stale manuscript/docs language from 25k pilot framing to scaled latent and mapped-external framing.
- [ ] Decide whether to retire or refresh the legacy `workflows/Snakefile` surface.
- [ ] Group future large-run outputs under explicit `data/processed/latent/`, `data/processed/external/`, and `results/tables/latent/` namespaces if/when a larger breaking path cleanup is acceptable.

## 1. Gateway scANVI/scArches Reference Mapping

- [x] Build a Gateway reference training script using existing scaled/lineage scVI inputs where possible.
- [x] Train scANVI labels from Gateway fine/coarse cell states.
- [x] Save reference model artifacts, label vocabulary, gene set, and registry metadata.
- [x] Implement GSE184117 query preparation with exact gene intersection and metadata harmonization.
- [x] Project GSE184117 into Gateway latent/reference space using scArches-style query mapping.
- [x] Emit external query latent coordinates, predicted labels, label confidence/entropy, donor features, and mapping QC.
- [x] Compare scANVI/scArches features with current marker-reference features.
- [x] Update report/manuscript only if mapping confidence and concordance improve or clarify limitations.

## 2. Full 4M-Cell Latent Model

- [x] Estimate local memory/runtime from the 250k run and decide local M5 versus remote server execution.
- [x] Add remote-run notes for the T1000/enterprise server path, including environment, data transfer, resumability, and output sync.
- [x] Add reproducible full-4M scVI Make targets and provenance entries.
- [x] Add SSH/rsync/tmux remote runner for `sabharwaly2@mia.ninds.nih.gov`.
- [x] Launch first all-cell 3,003-gene attempt on `mia` and diagnose silent preprocessing stop before GPU training.
- [x] Add a 1M-cell stratified fallback target after all-cell materialization proved memory-limited on `mia`.
- [x] Add a 500k-cell stratified target after the 1M sampled materialization also approached the server memory ceiling.
- [x] Add a second 250k-cell seed-stability target as the known-safe scale for the current H5AD/scVI materialization route.
- [x] Add a chunked all-cell reduced-H5AD builder so 4M-cell training can be attempted without the memory-hostile raw backed row/column slice.
- [x] Train full Gateway scVI/scANVI model or a defensible sketch/online alternative.
- [x] Validate full-model label purity, technical mixing, donor mixing, marker continuity, and seed stability.
- [ ] Compare full-model embeddings with 250k and lineage-focused embeddings.

## 3. Milo Differential Abundance

- [x] Add R/Python workflow that consumes validated latent coordinates and donor metadata.
- [x] Define age bins and continuous-age neighborhood tests.
- [x] Run first-pass Milo-style DA on healthy-donor latent neighborhoods.
- [x] Run Milo-style DA on basal-to-neuronal and supporting/secretory neighborhoods.
- [x] Run publication-scale Milo-style DA on the full 4M reduced scVI atlas.
- [ ] Add matched technical sensitivity and donor-yield sensitivity.
- [ ] Summarize neighborhoods by marker enrichment and ORA feature themes.

## 4. Pseudotime / Palantir / CellRank

- [ ] Use lineage-focused scVI/scANVI embeddings as the substrate.
- [ ] Run diffusion pseudotime or Palantir-style fate probabilities for basal-to-OSN and basal-to-secretory paths.
- [ ] Evaluate whether RNA velocity or transition kernels are defensible before using CellRank.
- [ ] Test age/ORA shifts along pseudotime and fate probabilities.
- [ ] Require seed stability before manuscript claims.

## 5. cNMF / Program Discovery

- [ ] Run cNMF within HBC/GBC/progenitor cells.
- [ ] Run cNMF within immature/mature OSN cells.
- [ ] Run cNMF within sustentacular/secretory cells.
- [ ] Run cNMF within immune/stromal compartments.
- [ ] Test program scores against age, ORA residuals, chemistry/device, donor yield, and GSE184117 mapped states.
- [ ] Promote only seed-stable and technically robust programs.

## 6. Expanded External Validation Stack

- [ ] Query GEO/SRA/CELLxGENE for human olfactory/nasal single-cell and spatial datasets.
- [ ] Add a candidate registry table with accession, tissue, disease/age context, assay, raw availability, labels, donors, and validation role.
- [ ] Add GSE151973 bulk deconvolution using Gateway-derived signatures.
- [ ] Search COVID/anosmia/nasal biopsy datasets for regeneration-module validation.
- [ ] Contact/resolve original labels for GSE184117 if public annotations remain insufficient.

## 7. Accuracy / Model Improvement

- [ ] Add elastic-net stability paths.
- [ ] Add Bayesian additive regression tree exploration or document blocker.
- [ ] Add TabPFN/AutoGluon exploratory benchmarks with nested donor-level CV only.
- [ ] Add stacked ensembles only under nested CV.
- [ ] Prototype multitask/adversarial age prediction with chemistry/device removal, labeled as exploratory.
- [ ] Compare accuracy gains against interpretability and external concordance.

## 8. Biological Mechanism Analyses

- [ ] Add ligand-receptor or NicheNet-ready exports for immune/stromal/sustentacular-to-basal signaling.
- [ ] Test whether inflammatory/supporting signals plausibly regulate basal-cell stress or blocked neurogenesis.
- [ ] Cross-check mechanism candidates against GSE184117 stem-cell inflammation findings.
- [ ] Promote mechanism claims only if robust to cell-state, donor, and technical sensitivity.

## 9. Spatial / Histology Validation

- [ ] Search for human olfactory spatial transcriptomics or histology validation datasets.
- [ ] If spatial data exist, test cell2location-like mapping of Gateway states into neuroepithelium versus respiratory metaplasia.
- [ ] Draft wet-lab validation panel: TP63/KRT5, OMP, GAP43/DCX, CYP2A13, MUC5B, PTPRC/LST1.
- [ ] Define expected direction for each marker from ORA and external validation.

## Running Notes

- 2026-06-16: Created active tracker after repository scan. Generated cache directories were removed. Stale docs/manuscript text was updated to reflect scaled 250k/100k scVI and GSE184117 mapped-feature validation. Task 1 is now the active implementation focus.
- 2026-06-16: Added the first scANVI/scArches mapping workflow, including reference training/query projection CLI, Make target, command-manifest entry, QC/donor-feature helper functions, and unit tests. Focused tests, full tests, and ruff are clean before launching the heavier reference-mapping run.
- 2026-06-16: Completed Task 1 first pass. Trained a Gateway scANVI reference from `gateway_scvi_stratified_250k.h5ad`, reused the saved model for GSE184117 query transfer, and wrote `gse184117_scanvi_mapped.h5ad` with 59,656 query cells x 3,003 genes. All six external samples were high-confidence by label confidence/entropy. Same-named scANVI donor-feature concordance against Gateway age associations produced 58 rows: 27 concordant, 21 discordant, and 10 not evaluable, so the external claim remains small-n and mixed rather than independently validated.
- 2026-06-16: Started Task 2 compute setup. Non-interactive SSH to `sabharwaly2@mia.ninds.nih.gov` reaches the NIH host but currently fails authentication, so remote execution is prepared but not launched. Added full-4M Make/provenance targets, model saving in the scVI runner, a remote SSH/rsync/tmux launcher, and `docs/full_4m_compute_plan.md`. Recommendation is remote-first once SSH auth is available; local M5 remains a fallback only for smaller sketch runs.
- 2026-06-17: Launched the first all-cell 3,003-gene scVI attempt manually on `mia` through tmux after password-based SSH setup. The process reached 76 GB RSS during backed H5AD materialization, exited before GPU training, and left only the seed line in the log. Added unbuffered progress logging, skipped redundant HVG recomputation for fixed gene-list runs, and registered a marker-preserving 1,500-gene all-cell fallback target for memory-constrained full-cell training.
- 2026-06-17: The marker-preserving 1,500-gene all-cell fallback reached 103 GB RSS while still materializing the backed H5AD slice, so it was stopped before exhausting the 125 GB no-swap server. Added `scvi-scaled-1m`, a stratified 1M-cell run with the full 3,003-gene feature set, as the next defensible large-scale latent atlas step.
- 2026-06-17: The 1M-cell stratified run reached 106 GB RSS during sampled backed-slice materialization before training, so it was also stopped to protect the no-swap server. Added `scvi-scaled-500k` as the next memory-safe scale step with the full 3,003-gene feature set. True all-cell training now requires a chunked/Zarr-backed path or a higher-memory host rather than the current in-memory AnnData materialization route.
- 2026-06-17: The 500k-cell stratified run also reached 106 GB RSS during materialization before training and was stopped. Added `scvi-scaled-250k-seed23` to prioritize seed stability at the known-safe 250k scale while a chunked/Zarr-backed or higher-memory all-cell strategy is designed. The seed-stability target intentionally avoids fixed backed-column slicing because that HDF5 CSR path caused high transient memory even at 250k cells.
- 2026-06-17: Added `scripts/build_reduced_h5ad.py` and Make targets `scvi-reduced-4m`/`scvi-full-4m-reduced`. The new path writes 3,003-gene chunks from all Gateway rows and concatenates them on disk before all-cell scVI training, preserving all cells while avoiding the failed raw-H5AD backed fancy-index materialization route.
- 2026-06-18: Completed the all-cell reduced Gateway scVI run on `mia`. The chunked reducer wrote `data/processed/gateway_hvg3003_4m.h5ad` with 4,028,275 cells x 3,003 genes, and scVI trained on all cells for 100 epochs, writing `data/processed/gateway_scvi_full_4m_reduced.h5ad` with finite 10-dimensional `X_scvi` for all 4,028,275 cells. Validation on a 100k-cell sample passed embedding, fine/coarse label-purity, FLEX/device/condition/sex mixing, and basal/immature OSN/mature OSN/sustentacular marker-continuity checks. Progenitor and immune marker-continuity were limited in the full model, while the second 250k seed model passed immune continuity; keep mechanism claims gated until targeted lineage/program checks compare these models.
- 2026-06-18: Added and ran a Python Milo-style latent-neighborhood DA pilot on the local 250k scVI atlas. The workflow tested 1,000 healthy-donor neighborhoods with donor-level logit-fraction regressions adjusted for sex, chemistry, and collection method. All 1,000 neighborhoods were testable, but 0 passed BH FDR < 0.10 for age. The strongest nominal neighborhood was mostly quiescent HBC (`p=0.00568`, `FDR=0.885`), so this supports keeping neighborhood biology exploratory until sensitivity runs, full-model transfer, or stronger neighborhood definitions improve signal.
- 2026-06-18: Added targeted Milo-style lineage and secretory neighborhood runs. The lineage-focused 100k scVI pass tested 1,000 basal/INP/iOSN/mOSN/sustentacular neighborhoods and found 0 age neighborhoods at BH FDR < 0.10; the strongest nominal neighborhoods were early mature mOSN, late iOSN, and cycling HBC with negative age coefficients. The secretory/sustentacular/glandular pass also found 0 FDR-significant neighborhoods; nominal top neighborhoods were serous, club, mucous gland, and proliferating secretory with negative age coefficients. Treat this as hypothesis-generating support for regenerative/secretory depletion, not a promoted Milo claim.
- 2026-06-18: Upgraded the Milo-style workflow for publication-scale 4M runs. New `milo-full-4m`, `milo-full-4m-lineage`, and `milo-full-4m-secretory` targets use the all-cell reduced scVI atlas, 20,000 stratified seed neighborhoods, 100 nearest cells per neighborhood, minimum 30 donors, categorical technical covariates, and numeric donor-yield adjustment. The 250k/100k targets are now explicitly engineering/sensitivity runs, not final scientific endpoints.
- 2026-06-18: Completed the publication-scale full 4M Milo-style runs on `mia` in roughly 17 minutes. The broad all-cell run tested 20,000 neighborhoods and found 3,673 with age FDR < 0.10, mostly negative age coefficients; the top neighborhood was mucous gland/respiratory secretory (`FDR=1.62e-05`). The lineage-focused full 4M run found 5,613 age neighborhoods at FDR < 0.10, dominated by negative early iOSN, late iOSN, INP, HBC/suprabasal, and mature OSN neighborhoods. The secretory-focused full 4M run found 285 neighborhoods at FDR < 0.10 but with mixed direction inside the restricted secretory denominator. Keep matched sensitivity and marker/program annotation as the next claim gate.
