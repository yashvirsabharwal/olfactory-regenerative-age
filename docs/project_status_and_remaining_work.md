# ORA Research Tracker

Updated: 2026-06-25

This is the single detailed tracker for the ORA research program. The focus is the science: what evidence exists, what claims are supported, what limitations remain, and what research tasks are left.

## Current Research Position

ORA currently supports a conservative but real research claim: healthy human olfactory epithelium contains a modest, reproducible, donor-level regenerative aging tissue-state axis. The project is strongest as an interpretable atlas-scale reanalysis of olfactory epithelial aging, not as an absolute age clock, disease diagnostic, lineage-flux measurement, or causal mechanism paper.

The remaining research gaps are mostly about strengthening validation and biological interpretation:

- stronger independent or orthogonal validation;
- perturbation/organoid adapter evidence;
- clinical/odor phenotype coupling;
- heavier external-prior regulon and ligand-receptor validation;
- release-grade atlas product definition.

## Supported Claims

| Claim | Status | Evidence | Boundary |
| --- | --- | --- | --- |
| Healthy-donor ORA aging axis | Supported | Repeated donor-level CV, shuffled-age null, calibration, sensitivity analysis, model card, feature interpretation. | Modest regenerative tissue-state axis, not an absolute age clock. |
| Composition-aware aging signal | Supported | Donor-level CLR compositional modeling and strict/single-context sensitivity. | Association only; not lineage flux. |
| Technical confounding addressed | Supported with caveat | Negative controls, technical-only/yield-only models, backend provenance, leave-context-out tests. | Technical covariates matter; model is not context-invariant. |
| External validation | Limited support | GSE184117 public-matrix reanalysis, scANVI/mapped features, public-data exhaustion, author response. | Small-n/contextual only; no author-label replication. |
| Cross-tissue specificity | Supported with caveat | CELLxGENE nasal/bronchial donor-level age-effect estimates plus first-pass ORA feature classification. | Directional/contextual; no FDR-significant adult comparator effects yet and LungMAP is developmental-context only. |
| Spatial/histology validation | Designed, data gap recorded | Public spatial triage, marker panel, readout plan, and targeted validation design. | No direct public adult human olfactory spatial aging dataset found; experimental validation needed. |
| Perturbation/organoid search | Designed, adapter candidates found | Public organoid/ALI perturbation triage and minimum experiment design. | No clean adult olfactory aging perturbation evidence yet; GSE309325/GSE299529 require adapters/audits. |
| Foundation/aging-clock comparison | Supported | scVI donor embedding, ORA+scVI hybrid, Geneformer V1, expression-clock baseline. | Best strict expression-clock MAE is 11.98 years, not below 10. |
| Regeneration-axis interpretation | Supported as annotation | Feature map, pathway modules, age associations, ORA correlations. | Pathway-prior association, not causality. |
| Regulatory drivers | Hypothesis layer | Curated TF/pathway target-program scoring. | Not SCENIC/decoupler regulon inference. |
| Niche ligand-receptor signaling | Hypothesis layer | Curated sender-receiver pseudobulk LR scoring. | Not spatial colocalization, receptor activation, or causal niche signaling. |
| Latent neighborhood remodeling | Supported with limitation | Full 4M scVI/Milo-style analysis, edgeR parity, age-bin checks, official MiloR subset sensitivity. | Remodeling/density only; not fate. |
| Regeneration topology | Exploratory | DPT pilot on basal-to-neuronal lineage subset. | No velocity layers; no lineage flux/rate claim. |
| AD/PD projection | Exploratory only | Frozen healthy-trained ORA projection and label permutation. | No diagnostic, biomarker, or disease-mechanism claim. |
| Disease DE | Hypothesis-generating | edgeR/limma/matched/sentinel audits. | Method-sensitive; not primary biology. |

## Research Completed

### 1. Code, Data, And Claim Cleanup

Status: complete.

- Corrected lineage-ratio names to match formulas.
- Hardened age parsing for mixed numeric/string metadata.
- Made cohort inclusion rules explicit.
- Added model backend provenance and fallback gating.
- Created claim boundaries and prohibited overclaims.

Key outputs:

- `data/processed/cohort_manifest.tsv`
- `data/processed/ora_feature_matrix.tsv`
- `data/processed/ora_augmented_feature_matrix.tsv`
- `docs/claim_ledger.md`

### 2. Core ORA Model

Status: complete.

- ORA is trained on healthy donors only.
- Repeated donor-level CV and shuffled-age null show signal above null.
- Calibration and sensitivity outputs exist.
- Feature interpretation maps model signal to cell-state, module, and lineage features.

Key outputs:

- `results/tables/ora_repeated_cv_summary.tsv`
- `results/tables/ora_permutation_empirical.tsv`
- `results/tables/ora_calibration.tsv`
- `results/tables/ora_feature_interpretation.tsv`
- `results/tables/ora_model_card.tsv`

### 3. Statistical Hardening

Status: complete.

- Donor-level CLR compositional age models are complete.
- Negative-control models show technical covariates are important but do not fully explain the best ORA signal.
- Feature-family ablation identifies expression PCs and hybrid ORA/scVI features as strong comparators.
- Leave-one-context-out tests show ORA is not fully context-invariant.

Key outputs:

- `results/tables/compositional_age_model_summary.tsv`
- `results/tables/compositional_age_model_sensitivity.tsv`
- `results/tables/ora_negative_control_performance.tsv`
- `results/tables/feature_family_ablation_model_comparison.tsv`
- `results/tables/leave_context_out_summary.tsv`

### 4. External Validation And Public Data Exhaustion

Status: complete with limitation.

- GSE184117 is the only direct public human olfactory aging/presbyosmia scRNA-seq validation target found.
- Usable external comparison remains 3 healthy versus 3 presbyosmic biopsy donors.
- Author response received 2026-06-25 confirmed:
  - raw sequencing files and 10x matrices are public through GEO;
  - sample names and clinical metadata can be matched through GEO records and the study's Table 1;
  - no separate per-cell manual annotation table exists.
- Therefore, GSE184117 is public-matrix reanalysis/contextual support, not author-label replication.

Key outputs:

- `results/reports/gse184117_reanalysis_status.md`
- `results/reports/public_data_exhaustion.md`
- `results/tables/external_candidate_matrix.tsv`
- `results/tables/public_data_exhaustion.tsv`
- `results/tables/external_validation_evidence.tsv`
- `docs/external_label_request_log.md`

### 5. Foundation And Aging-Clock Benchmarks

Status: complete.

- Fixed Gateway benchmark subsets were generated.
- scVI donor embedding benchmark is complete.
- ORA+scVI hybrid benchmark is complete.
- Geneformer V1-10M benchmark is complete and weaker than ORA/scVI/expression baselines.
- Fold-internal expression-clock baseline is complete. Best strict expression-clock model reaches MAE 11.98 years.

Research conclusion:

- ORA is not the strongest pure age-prediction model, but it is more interpretable and biologically structured.
- Current evidence does not justify a below-10-year MAE claim.

Key outputs:

- `docs/foundation_model_benchmark_plan.md`
- `results/tables/foundation_model_benchmark.tsv`
- `results/tables/foundation_model_runtime.tsv`
- `results/tables/geneformer_age_model_summary.tsv`
- `results/tables/aging_clock_baseline_model_comparison.tsv`

### 6. Cross-Tissue Specificity

Status: complete with caveat.

Completed:

- Built a cross-tissue comparator plan and candidate matrix.
- Classified ORA/module features as olfactory-specific, airway/nasal shared, pan-epithelial regenerative, immune/inflammatory shared, or not comparable.
- Downloaded and analyzed three public CELLxGENE comparator H5ADs:
  - Xu nasal scRNA-seq: 34,833 cells, 7 adult donors, ages 30-66.
  - Xu bronchial scRNA-seq: 2,075 cells, 17 adult donors, ages 51-87.
  - LungMAP human lung age-groups: 46,500 cells, 9 donors/stages, but only 3 adult donors all at age 31.
- Computed donor-level comparator composition/module age effects.
- Generated 71 adult nasal/bronchial comparator estimates.
- Mapped 74 ORA features to measured adult comparator age-effect rows.

Main result:

- Nasal/bronchial comparators now support explicit shared-airway and shared-immune context for many ORA features.
- Adult comparator estimates are mostly directional: 3 nominal p<0.05 effects and 0 within-dataset/scope FDR<0.05 effects.
- LungMAP is useful only as developmental context unless more adult lung donors are added.
- Non-airway specificity against skin/gut/blood remains desirable, but the airway/nasal specificity blocker is now cleared.

Key outputs:

- `docs/cross_tissue_specificity_plan.md`
- `docs/cross_tissue_age_effects.md`
- `results/tables/cross_tissue_candidate_matrix.tsv`
- `results/tables/cross_tissue_cellxgene_asset_inventory.tsv`
- `data/processed/cross_tissue_cellxgene_donor_features.tsv`
- `results/tables/cross_tissue_cellxgene_module_coverage.tsv`
- `results/tables/cross_tissue_age_effects.tsv`
- `results/tables/ora_cross_tissue_age_effect_summary.tsv`
- `results/tables/ora_cross_tissue_feature_classification.tsv`
- `results/tables/ora_cross_tissue_specificity.tsv`

### 7. Regeneration-Axis Biology

Status: complete as interpretation layer.

- Created controlled regeneration-axis feature map.
- Added curated regeneration pathway modules.
- Tested module age associations and ORA correlations.

Main biological pattern:

- Immature-OSN and mature-OSN pathway scores decline significantly with age in primary adjusted models.
- Strict-cohort models show broader context-sensitive negative shifts across neurogenic, repair, inflammatory, oxidative-stress, senescence, and respiratory-metaplasia programs.

Key outputs:

- `resources/feature_maps/regeneration_axis_feature_map.tsv`
- `results/tables/regeneration_axis_feature_map.tsv`
- `results/tables/regeneration_axis_theme_summary.tsv`
- `configs/regeneration_gene_sets.yaml`
- `results/tables/regeneration_module_age_associations.tsv`
- `results/tables/regeneration_module_ora_correlations.tsv`

### 8. Regulatory Drivers

Status: complete as hypothesis layer.

- Curated TF/pathway target-program scoring was run on genomewide pseudobulk.
- Top ranked driver hypotheses include WNT/CTNNB1, ASCL1, STAT/IRF, TP63, OMP/ADCY3, NOTCH/HES, YAP/TEAD, and NF-kB/IL17/TNF.

Boundary:

- These are curated target-program hypotheses, not inferred regulons and not causal evidence.

Key outputs:

- `configs/regulatory_driver_gene_sets.yaml`
- `docs/regulatory_driver_feasibility.md`
- `results/tables/regulatory_driver_map.tsv`
- `results/tables/regulatory_driver_age_associations.tsv`
- `results/tables/regulatory_driver_ora_correlations.tsv`

### 9. Niche Ligand-Receptor Signaling

Status: complete as hypothesis layer.

- Curated ligand-receptor pseudobulk scoring was run by donor, sender group, and receiver group.
- Twelve families were tested: IL17, TNF, IFN, IL6, MIF, CXCL, CCL, EGFR/AREG, Notch, Wnt, TGF-beta, and IL1.
- Gene coverage is strong: minimum 88.9%, mean 98.4%.
- Analysis produced 124 sender-receiver hypotheses and 248 age-association rows.

Main pattern:

- Top hypotheses include immune-to-lineage IFN, TNF, CCL/CXCL, and MIF signals.
- Many inflammatory-prior signals decline with age/ORAA, so the current story is context-sensitive niche remodeling rather than simple increased inflammaging.

Boundary:

- These are pseudobulk sender-receiver hypotheses, not physical cell-cell communication, spatial colocalization, receptor activation, or causal niche control.

Key outputs:

- `configs/niche_ligand_receptor_pairs.yaml`
- `docs/niche_signaling_feasibility.md`
- `results/tables/niche_ligand_receptor_gene_coverage.tsv`
- `results/tables/niche_ligand_receptor_donor_scores.tsv`
- `results/tables/niche_ligand_receptor_age_associations.tsv`
- `results/tables/niche_ligand_receptor_ora_associations.tsv`
- `results/tables/niche_driver_priority_table.tsv`

### 10. Latent Neighborhoods

Status: complete with limitation.

- Full 4M scVI/Milo-style lineage-neighborhood analyses are complete.
- Matched and all-donor results are available.
- EdgeR parity, age-bin robustness, program scoring, and official MiloR subset sensitivity are complete.

Main pattern:

- Age-associated lineage-neighborhood remodeling is supported.
- Matched analysis narrows the strongest exact-neighborhood signal to Early iOSN/iOSN-related structure.

Boundary:

- This supports neighborhood remodeling/density shifts, not fate, trajectory, or lineage flux.

Key outputs:

- `results/tables/scvi_embedding_claim_gates.tsv`
- `results/tables/milo_full_4m_lineage_summary.tsv`
- `results/tables/milo_full_4m_lineage_matched_summary.tsv`
- `results/tables/milo_full_4m_lineage_edger_parity_summary.tsv`
- `results/tables/milor_lineage_subset_summary.tsv`
- `results/tables/milo_full_4m_lineage_matched_program_summary.tsv`

### 11. Fate/Regeneration Dynamics

Status: complete with limitation.

- Input audit found no `spliced`/`unspliced` layers in inspected H5ADs.
- RNA velocity and CellRank velocity-kernel fate inference are no-go.
- Scanpy DPT pilot on 40,000 basal-to-neuronal lineage cells found exploratory positive ordering between expected lineage rank and DPT pseudotime.

Boundary:

- DPT supports topology only. It does not support lineage flux, fate transition rate, or causal regeneration dynamics.

Key outputs:

- `results/tables/regeneration_dynamics_input_audit.tsv`
- `results/tables/regeneration_dynamics_feasibility.tsv`
- `results/reports/regeneration_dynamics_feasibility.md`
- `results/tables/regeneration_dynamics_summary.tsv`

### 12. Disease Projection And Disease DE

Status: exploratory/deferred.

- AD/PD ORA projections exist but disease groups are small and confounded by technical context.
- Disease DE audits exist but remain method-sensitive.

Boundary:

- No AD/PD diagnostic, biomarker, or disease-mechanism claim.
- Disease observations are exploratory stress tests only.

Key outputs:

- `results/tables/ndd_ora_projection_summary.tsv`
- `results/tables/ndd_label_permutation.tsv`
- NDD guardrail and DE audit summary tables in `results/tables/`

### 13. Spatial And Histology Validation Design

Status: complete as design; public data gap remains.

Completed:

- Re-checked public human olfactory/nasal spatial resources.
- Confirmed no direct public adult human olfactory epithelial spatial aging dataset was found.
- Triage table records:
  - GSE235714 as human nasal CRS/healthy-control GeoMx spatial context only.
  - GSE292993 as airway/lung spatial context only.
  - GSE303809 as fetal olfactory/head-section MERFISH developmental context only.
- Defined a 9-panel marker set covering HBC, activated HBC/INP, iOSN, mOSN, sustentacular/support, respiratory metaplasia, immune/inflammatory, stress/senescence, and niche LR signals.
- Defined donor-level spatial readouts and a target adult olfactory histology/spatial experiment.

Boundary:

- This is a validation design, not measured orthogonal evidence yet.
- It can support future localization claims but cannot prove lineage flux, regeneration rate, or causality.

Key outputs:

- `docs/spatial_perturbation_validation_plan.md`
- `configs/spatial_validation_markers.yaml`
- `results/tables/spatial_validation_candidate_matrix.tsv`
- `results/tables/spatial_validation_marker_panel.tsv`
- `results/tables/spatial_validation_readout_plan.tsv`
- `results/tables/spatial_validation_search_log.tsv`

### 14. Perturbation, Organoid, And ALI Search/Design

Status: complete as search/design; measured adapter evidence still needed.

Completed:

- Searched public GEO/web resources for human olfactory organoid, olfactory injury/regeneration, nasal epithelial cytokine, IFN/TNF, and ALI perturbation datasets.
- Found no clean adult human olfactory aging perturbation dataset.
- Ranked seven public candidates:
  - GSE309325: highest-priority olfactory-relevant human nasal organoid scRNA-seq infection/time-course lead.
  - GSE299529: strongest nasal ALI TNF/TGF-beta scRNA-seq perturbation lead.
  - GSE286616 and GSE309353: bulk nasal ALI viral/IFN/NF-kB context candidates.
  - GSE175541, GSE324335, and GSE271245: lower-priority context-only candidates.
- Defined a minimum adult olfactory organoid/ALI perturbation experiment.

Boundary:

- This is not yet measured perturbation support.
- Current public candidates can strengthen plausibility only after adapters score ORA/regeneration modules by condition/timepoint.

Key outputs:

- `docs/perturbation_validation_plan.md`
- `configs/perturbation_validation_candidates.yaml`
- `results/tables/perturbation_validation_candidate_matrix.tsv`
- `results/tables/perturbation_validation_search_log.tsv`
- `results/tables/perturbation_minimum_experiment_design.tsv`

## Research Still Left

### Priority 1: Perturbation Adapter Implementation

Why it matters:

- This is the best route from association to causal support.

Tasks:

- Audit and attempt a GSE309325 organoid adapter first.
- If GSE309325 is blocked, try GSE299529 Seurat RDS conversion and module scoring.
- If single-cell adapters are blocked, implement a bulk module adapter for GSE286616 or GSE309353.
- Score ORA/regeneration/metaplasia/IFN/TNF/IL17/senescence modules by condition and timepoint where possible.
- Keep all claims guarded as perturbation-context evidence unless adult olfactory donor-level responses reproduce ORA directions.

Candidate public adapters:

- GSE309325: human nasal organoids with olfactory and respiratory epithelium, mock versus SARS-CoV-2 time course.
- GSE299529: nasal ALI TNF/TGF-beta scRNA-seq, processed Seurat RDS.
- GSE286616: nasal ALI rhinovirus plus IRF3/NF-kB inhibitor bulk RNA-seq.
- GSE309353: primary nasal ALI RSV bulk RNA-seq.

Done when:

- At least one perturbation adapter emits module-level contrast tables, or all top candidates have documented hard blockers.

### Priority 2: Clinical/Odor Phenotype Coupling

Why it matters:

- Chronological age alone is biologically limited. Odor function, smoking, biopsy site, CRS/COVID history, or UPSIT/SIT would make ORA more clinically meaningful.

Tasks:

- Audit available Gateway metadata again for smell/clinical fields.
- Use GSE184117 sample metadata/SIT values only as small-n context.
- Identify whether any public olfactory dataset has usable odor scores with expression.
- Keep this blocked unless real covariates are available.

Done when:

- `results/tables/clinical_coupling_summary.tsv` exists, or the project records a clear blocked/no-data decision.

### Priority 3: Heavier External-Prior Mechanism Validation

Why it matters:

- Current regulatory-driver and LR layers are transparent curated fallbacks. Heavier methods would make mechanism claims more SOTA.

Tasks:

- Try decoupler with CollecTRI/DoRothEA or PROGENy for pathway activity.
- Consider pySCENIC/SCENIC+ only if runtime and inputs are practical.
- Try LIANA/NicheNet/CellPhoneDB/CellChat only if package/database setup is stable.
- Compare heavy-method results to current curated drivers/LR hypotheses.

Done when:

- `results/tables/external_prior_regulatory_driver_validation.tsv` or equivalent exists, or a documented no-go exists.

### Priority 4: Atlas Product Definition

Why it matters:

- A release-grade atlas object would make this a stronger reusable research resource.

Tasks:

- Define what can be redistributed.
- Define release H5AD schema.
- Include ORA scores, model features, module scores, scVI/foundation embeddings, validation mappings, and metadata dictionary where allowed.
- Record checksums and source provenance.

Done when:

- `docs/atlas_release_schema.md` and a metadata dictionary exist.

### Priority 5: New External Cohorts

Why it matters:

- This is the only clean way to convert external validation from limited/contextual to strong.

Tasks:

- Continue watching GEO/SRA/CELLxGENE for human olfactory aging, presbyosmia, anosmia, CRS, COVID smell-loss, and AD/PD olfactory biopsy datasets.
- Promote only datasets with donor-level age and suitable olfactory epithelial expression.
- Keep GSE184117 as the best current direct public validation but do not overstate it.

Done when:

- A stronger external dataset is found and processed, or the tracker records the current no-stronger-public-data state.

## Research Priority Order

1. Perturbation/organoid/ALI search and design.
2. Clinical/odor phenotype coupling audit.
3. Heavier external-prior regulatory/LR validation.
4. Atlas release schema.
5. New external cohort monitoring.

## Current Verification Snapshot

Latest clean verification:

- `ruff check .`: passed.
- `make cross-tissue-age-effects`: passed; 246 effect rows, 71 adult nasal/bronchial estimates.
- `make spatial-validation-plan`: passed; 4 candidate rows, 9 marker panels, 5 readout families.
- `make perturbation-validation-plan`: passed; 7 candidates, 2 high-priority adapter leads.
- `make test`: passed, 141 tests.
- `make output-provenance`: passed, 111 stages, 316 outputs, 0 missing.

## Bottom Line

The core research is now strong for an interpretable atlas-scale study of olfactory regenerative aging. What remains scientifically is not more basic model hardening; it is orthogonal validation beyond public transcriptomic reanalysis: perturbation adapter evidence, clinical phenotype coupling, heavier mechanism-method validation, and a release-grade atlas product.
