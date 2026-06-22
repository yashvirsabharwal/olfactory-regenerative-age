# ORA Publication Readiness Tracker

Updated: 2026-06-22

This is the active task board for moving ORA from a strong internal analysis to a defensible preprint/manuscript. It intentionally prioritizes claim gates, figure/manuscript readiness, and reproducibility over open-ended benchmarking.

## Publication Frame

Current target: a computational single-cell atlas reanalysis/resource paper centered on a healthy human olfactory regenerative aging axis.

Primary claim:

- Healthy donor olfactory epithelial composition and curated module features encode a modest, reproducible regenerative aging axis.

Secondary mechanistic claim:

- Full 4M scVI neighborhoods support a secondary layer of age-associated lineage-neighborhood remodeling. Exact Python-neighborhood analyses support a narrow matched Early iOSN/immature-neuron depletion result, enriched for immature-neuron genes and depleted for HBC programs; official MiloR subset sensitivity instead emphasizes broader HBC/suprabasal/sustentacular structure, so Early iOSN wording should remain guarded.

Exploratory only:

- GSE184117 concordance, AD/PD projection, genome-wide disease DE, secretory-only neighborhood DA, broad all-donor neighborhood maps, and any future pseudotime/cNMF/ligand-receptor analyses.

Do not claim:

- Absolute biological-age clock accuracy.
- AD/PD diagnostic utility.
- Measured lineage flux from cross-sectional data.
- UMAP-derived trajectory or neighborhood biology.

## Current Status Snapshot

- Repository is maintained on `main`; generated data/results remain intentionally ignored.
- Full Gateway analysis covers 202 donors and 4,028,275 cells.
- ORA modeling uses 187 healthy age-usable donors with donor-level CV and shuffled-age null testing.
- Full all-cell reduced scVI was trained on all 4,028,275 cells with 3,003 HVG/marker genes and finite 10-dimensional `X_scvi`.
- Full 4M Milo-style lineage result: 5,613 / 20,000 neighborhoods at age FDR < 0.10, dominated by negative neuronal/regenerative neighborhoods.
- Matched FLEX v2/device lineage result: 1 / 20,000 neighborhood at age FDR < 0.10, Early iOSN, `age_coef=-1.014`, `FDR=0.0427`.
- Matched Early iOSN neighborhood program scores: immature neuron `z=2.91`, senescence/SASP `z=1.63`, HBC activation/injury `z=-0.84`, HBC identity `z=-0.88`.
- Full 4M lineage age-bin robustness is complete: 4,705 / 5,613 negative age-associated lineage neighborhoods are lower in the oldest observed donor bin than the youngest observed donor bin. The strict matched Early iOSN hit also agrees directionally across observed bins (`lt45` to `60_74`; old-minus-young median logit fraction `-0.944`).
- edgeR count-model parity supports the Python Milo-style direction: signed-effect Spearman is `0.916` for all-donor lineage and `0.924` for matched lineage; the single matched Python-significant Early iOSN neighborhood overlaps the matched edgeR-significant set.
- Official MiloR subset sensitivity is complete and narrows the interpretation: it confirms strong age-associated latent-neighborhood structure, but dominant official MiloR signals are HBC/suprabasal/sustentacular rather than a direct matched Early iOSN replication.
- Run hierarchy is standardized in `docs/run_hierarchy.md`: full 4M reduced scVI is the primary latent substrate; full 4M Python Milo-style neighborhoods are the primary neighborhood map; edgeR parity and official MiloR subsets are sensitivity layers; 25k/100k/250k/500k/1M runs are pilot, reference, or optional stress-test runs.
- GSE184117 has scANVI/scArches mapping and donor-feature concordance, but remains small-n and mixed.
- NDD projection remains exploratory: 5 AD and 5 PD donors, all FLEX v2/device.
- Genome-wide DE has edgeR/limma parity and audits, but disease biology remains hypothesis-generating.

## Priority 1: Manuscript Claim Gates

These are the tasks that most directly determine whether the current story is ready to preprint.

- [x] Add age-bin robustness for full 4M lineage Milo-style neighborhoods.
  - Goal: show whether Early iOSN/immature-neuron depletion is monotonic or restricted to a specific age range.
  - Outputs: `results/tables/milo_full_4m_lineage_age_bin_neighborhoods.tsv`, `results/tables/milo_full_4m_lineage_age_bin_summary.tsv`, `results/tables/milo_full_4m_lineage_matched_age_bin_neighborhoods.tsv`, and `results/tables/milo_full_4m_lineage_matched_age_bin_summary.tsv`.
  - Result: broad all-donor negative lineage neighborhoods largely agree with age-bin direction; the single matched Early iOSN hit remains directionally negative from the youngest to oldest observed bins, though per-neighborhood bin donor counts are sparse.

- [x] Add Milo implementation-parity decision.
  - Goal: either run a focused official MiloR parity subset or write a transparent rationale for the Python donor-level neighborhood model.
  - Outputs: `docs/milo_implementation_parity.md`, exact-neighborhood edgeR parity tables, official MiloR full-lineage subset tables, and official MiloR matched-lineage subset tables.
  - Result: edgeR supports the exact Python-neighborhood age-effect directions, including the matched Early iOSN hit. Official MiloR confirms broad age-associated latent-neighborhood structure but does not make the matched Early iOSN result its dominant independent finding, so manuscript language should remain conservative.

- [ ] Compare full 4M scVI embedding with 250k seed and lineage-focused embeddings.
  - Goal: verify that the key lineage/matched signals are not an artifact of the full reduced model only.
  - Required outputs: embedding comparison summary, overlapping marker-continuity and theme-level concordance.
  - Acceptance: full 4M model is manuscript-primary; smaller models are clearly sensitivity anchors.

- [ ] Refresh claim ledger after the embedding comparison gate.
  - Goal: convert current exploratory/mechanistic language into final preprint language.
  - Required outputs: updated `docs/claim_ledger.md`, `docs/methodology_standards.md`, and manuscript result text.
  - Acceptance: every promoted claim has a matching audit/sensitivity table.

## Priority 2: Manuscript And Figure Package

These tasks turn the science into a coherent paper.

- [ ] Update manuscript framework for the new full 4M neighborhood and program-enrichment results.
  - Add matched Early iOSN program enrichment to the Results structure.
  - Reframe full 4M neighborhood analysis as a mechanistic support layer, not a central standalone paper.

- [ ] Update the LaTeX manuscript draft.
  - Add current methods for full 4M scVI, Milo-style DA, matched sensitivity, membership export, and program scoring.
  - Update Results, Discussion, Limitations, and figure callouts.
  - Keep citations verified and DOI-backed.

- [ ] Redesign main and extended figures.
  - Figure 1: cohort and claim-gated workflow.
  - Figure 2: ORA composition/module aging axis.
  - Figure 3: model performance, null, calibration.
  - Figure 4: stable ORA feature biology.
  - Figure 5: external validation and NDD guardrails.
  - Figure 6: full 4M scVI/Milo/program enrichment.
  - Extended Data: audits, model card, external evidence, DE parity, scVI validation, full neighborhood tables.

- [ ] Generate publication tables.
  - Required: cohort summary, model card, top stable features, external evidence ledger, NDD appendix, DE audit summary, full 4M Milo summary, matched program enrichment summary.

- [ ] Compile manuscript PDF once TeX tooling is available.
  - Acceptance: PDF builds from clean repo state and all figure/table paths resolve.

## Priority 3: External Validation Upgrade

This remains the largest scientific weakness.

- [ ] Run a fresh external dataset search for human olfactory/nasal single-cell, bulk, spatial, COVID/anosmia, presbyosmia, and aging datasets.
  - Required outputs: candidate registry with accession, tissue, assay, donors, age/disease metadata, raw availability, and validation role.
  - Acceptance: each dataset is labeled as direct validation, marker/context only, blocked, or not useful.

- [ ] Decide whether any dataset can provide donor-level ORA validation beyond GSE184117.
  - Direct validation requires age/status metadata, expression, and cell/sample labels or usable reference mapping.
  - If none exist, document this as a limitation rather than continuing open-ended searches.

- [ ] Add GSE151973 bulk/deconvolution context only if it clarifies olfactory vs respiratory marker specificity.
  - Acceptance: no claim that it validates donor-level ORA aging unless donor-level age design supports it.

- [ ] Consider author-label recovery for GSE184117.
  - Action: draft a concise data-request note for original cell labels or annotation files.
  - Acceptance: either labels obtained and integrated, or request documented as an unresolved limitation.

## Priority 4: Optional Mechanistic Extensions

These are valuable only after Priority 1 gates are done.

- [ ] Pseudotime or diffusion/Palantir-style lineage analysis.
  - Use only validated non-UMAP latent space.
  - Promote only if age/ORA shifts are stable and not chemistry/device driven.

- [ ] cNMF/program discovery in focused compartments.
  - Candidate compartments: HBC/GBC/progenitor, iOSN/mOSN, sustentacular/secretory, immune/stromal.
  - Promote only if programs are stable across seeds and interpretable beyond curated gene-set scoring.

- [ ] Ligand-receptor or NicheNet-ready export.
  - Focus on immune/stromal/sustentacular signals to basal or immature neuronal states.
  - Treat as hypothesis-generation unless externally supported.

- [ ] Spatial/histology validation search.
  - Search for human olfactory spatial transcriptomics or histology datasets.
  - If no public data exist, keep wet-lab validation panel as a proposed follow-up.

## Priority 5: Reproducibility And Repository Hardening

- [ ] Decide whether to retire or refresh `workflows/Snakefile`.
  - Acceptance: no stale workflow suggests an outdated MVP-only path.

- [ ] Add a manuscript rerun profile.
  - Could be a Make aggregate target, Snakemake profile, or documented staged rerun instructions.
  - Acceptance: user can regenerate manuscript-facing outputs without guessing command order.

- [ ] Preserve remote-compute notes.
  - Include `mia` tmux/rsync workflow, memory notes, runtime notes, and which large artifacts stay remote.

- [ ] Update output provenance after final figure/table generation.
  - Acceptance: 0 missing outputs for manuscript-facing artifacts.

- [ ] Run final checks before preprint.
  - `.venv/bin/python -m ruff check .`
  - `PYTHON=.venv/bin/python make test`
  - `PYTHON=.venv/bin/python make output-provenance`
  - citation-key check
  - figure existence check
  - PDF build check

## Completed Major Milestones

- [x] Gateway H5AD inspection and cohort manifest.
- [x] Donor-level composition, CLR, lineage-ratio, and module feature matrices.
- [x] Healthy-only ORA models with repeated donor-level CV.
- [x] Shuffled-age permutation/null testing.
- [x] Nested tuning and stacking benchmarks.
- [x] ORA calibration, residual diagnostics, sensitivity tables, and model card.
- [x] GSE184117 raw 10x module scoring, marker panels, mapped donor features, and scANVI/scArches query mapping.
- [x] AD/PD frozen ORA projection with uncertainty, matched context, and label permutation.
- [x] Targeted and genome-wide pseudobulk export.
- [x] Genome-wide edgeR and limma-voom summaries with matched and sentinel-gene audits.
- [x] Latent-space audit showing CELLxGENE exposes only UMAP.
- [x] Stratified 250k scVI, second seed target, lineage-focused scVI, and full 4M reduced scVI.
- [x] Full 4M scVI validation.
- [x] Full 4M Milo-style all-cell, lineage, secretory, matched, and donor-yield-adjusted analyses.
- [x] Milo-style theme annotation and claim gates.
- [x] Neighborhood membership export.
- [x] Full 4M lineage and matched-lineage curated program enrichment.
- [x] Methodology standards and claim ledger.
- [x] Initial manuscript framework and LaTeX draft.

## Near-Term Execution Order

1. Full 4M vs 250k/lineage embedding comparison.
2. Update claim ledger and LaTeX draft with the standardized run hierarchy.
3. Refresh figures, especially the new full 4M scVI/Milo/program figure.
4. External validation search refresh and final validation-strength table.
5. Final reproducibility appendix and output provenance.
6. PDF build and preprint package review.
