# ORA Claim Ledger

Updated: 2026-06-25

## Supported Primary Claim

- Healthy donor olfactory epithelial composition contains a reproducible, modest age-associated regenerative-state signal above shuffled-age null models.
- Donor-level CLR compositional age modeling supports the primary ORA age-association directions for the significant stable cell states and records strict-threshold plus single chemistry/collection sensitivities.
- Negative-control analyses show technical covariates are non-negligible: simple technical-only ridge controls are competitive with or better than simple biological ridge controls, but remain worse than the best repeated boosted ORA candidate and worse than the stratified age-shuffle null.
- Stable ORA features can be organized into supporting/secretory epithelial, immune/inflammatory, stress/senescence, neuronal-lineage, and regenerative/progenitor themes, with each row treated as associational rather than causal.
- The full ORA/module feature set can now be organized on a controlled regeneration-axis vocabulary with literature-anchored expected age directions, confidence labels, Gateway age associations, repeated-CV importance, and cross-tissue specificity status.
- Curated regeneration pathway modules support a guarded mechanistic interpretation: in primary covariate-adjusted donor models, immature-OSN and mature-OSN pathway scores decline significantly with age, while strict-cohort results show broader context-sensitive negative shifts across neurogenic, repair, inflammatory, oxidative-stress, senescence, and respiratory-metaplasia modules.
- First-pass regulatory-driver analysis nominates WNT/CTNNB1, ASCL1, STAT/IRF, TP63, OMP/ADCY3, NOTCH/HES, YAP/TEAD, and NF-kB/IL17/TNF as ranked hypotheses from curated pseudobulk target programs. This is driver-prior evidence, not inferred regulon activity or causality.
- Curated ligand-receptor pseudobulk scoring nominates immune/support/epithelial niche hypotheses, especially IFN, TNF, CCL/CXCL, MIF, Wnt, and TGF-beta sender-receiver pairs linked to basal/progenitor or OSN receivers. This is associational sender-receiver evidence; many inflammatory-prior signals decline with age/ORAA, so they should be described as context-sensitive niche hypotheses.

## Allowed Claim Language

- ORA is a modest donor-level regenerative tissue-state axis for healthy human olfactory epithelium.
- Formal donor-level compositional modeling supports ORA as a composition-aware association, not as a causal lineage-flux measurement.
- ORA models are trained on healthy donors only; AD and PD donors are reserved for frozen-model projection.
- ORA should be interpreted as a relative, under-dispersed tissue-state readout rather than an absolute chronological-age predictor.
- ORA has measurable technical dependence; claim language should emphasize donor-level biological features plus explicit technical sensitivity checks, not technical-covariate-free prediction.
- AD/PD ORA projections are exploratory stress tests of healthy-trained models, not disease classifiers or biomarkers.
- GSE184117 is small-n supportive/contextual evidence for reference-transfer and feature-direction checks. The study team confirmed on 2026-06-25 that no separate per-cell manual annotation table is available; sample metadata can be matched from GEO records plus manuscript Table 1, so use independent reanalysis/annotation rather than author-label replication.
- Public-data exhaustion as of 2026-06-24 found no stronger public independent donor-level human olfactory epithelial healthy-aging sc/snRNA validation dataset. GSE302937 and post-COVID/long-COVID olfactory datasets can be used as disease/smell-loss context, not healthy-aging replication.
- Full 4M scVI and Milo-style neighborhood analyses support guarded density/remodeling interpretations, not direct fate, pseudotime, or lineage-flux measurements.
- Matched feature-family ablation shows donor-level pseudobulk expression PCs are the strongest current predictive family, followed by ORA+scVI hybrid and scVI donor embeddings. Strict fold-internal pseudobulk-PC validation lowers this to an aging-clock-style ridge MAE of 11.98 years, so expression-derived age signal is real and competitive but still not a <10-year or general biological-clock claim.
- Geneformer V1-10M has now been benchmarked locally as a modern foundation-model comparator on the fixed lineage subset, but its donor-level age model is weaker (best MAE 15.04 years) than scVI-only, ORA+scVI, and the fold-internal expression-clock baseline. This supports benchmarking completeness, not a new performance claim.
- First-pass cross-tissue specificity classification separates ORA/module features into olfactory-specific, airway/nasal shared, pan-epithelial regenerative, immune/inflammatory shared, and not-comparable classes. This is claim-triage evidence only until selected comparator tissues have measured age-effect estimates.
- Regeneration-axis feature-map summaries can be used to describe which ORA features are basal/quiescent, basal-activated, cycling, progenitor, immature-OSN, mature-OSN, sustentacular/barrier, secretory, respiratory-metaplastic, immune/inflammatory, or stress/senescent, provided the text remains associational and not causal.
- Regeneration module results can connect ORA to specific OE pathways such as TP63/HBC, ASCL1/NEUROG1/NEUROD1 neurogenesis, GAP43/DCX immature OSN, OMP/ADCY3 mature OSN, Notch/Wnt/YAP/EGFR repair, interferon, senescence/SASP, oxidative stress, and respiratory metaplasia, but only as cross-sectional association and pathway-prior evidence.
- Regulatory-driver results can be described as ranked mechanistic hypotheses from curated TF/pathway target programs scored in donor-cell-state pseudobulk. Acceptable wording: "nominates", "prioritizes", or "is consistent with"; avoid "proves", "drives", or "causes".
- Niche ligand-receptor results can be described as curated sender-receiver hypotheses from donor-cell-state pseudobulk. Acceptable wording: "nominates", "prioritizes", or "is consistent with"; avoid "physical communication", "spatial colocalization", "receptor activation", or "causal niche control".
- Leave-one-context-out robustness gives a direct conservative answer to confounding concerns: site cannot be evaluated because it is single-level/missing, FLEX v2 holdout is feasible and performs well, but collection-method, sex, and several ethnicity holdouts show substantial error inflation. Do not describe ORA or hybrid models as context-invariant.

## Supported Engineering Claim

- The repository can generate donor-level features, ORA model benchmarks, module summaries, NDD projections, targeted pseudobulk DE, genome-wide pseudobulk export, edgeR and limma-voom summaries, and an MVP from ignored local artifacts.
- Lineage-ratio feature names are numerator/denominator explicit after the 2026-06-24 cleanup: `ratio__inp_to_activated_hbc` preserves the INP-over-activated-HBC formula, and `ratio__mature_mosn_to_iosn` preserves the mature-mOSN-over-iOSN formula.
- The primary ORA training rule is explicitly manifest-recorded as healthy donor plus valid age; strict total-cell, lineage-cell, and mature-neuron threshold flags are recorded for sensitivity analyses rather than silently changing the primary cohort.
- The foundation-model benchmark path is runnable locally for Geneformer V1-10M and records gene coverage, runtime, peak memory, checkpoint, donor feature outputs, and transparent scGPT/scFoundation deferred rows.
- The cross-tissue specificity path is runnable locally and emits a comparator candidate matrix, feature-level classification table, summary table, plan document, and figure, with explicit pending status for external comparator age effects.
- The regeneration-axis feature-map path is runnable locally and emits a static resource map, result-level evidence join, theme summary, and figure for all 120 ORA/module features.
- The regeneration-module path is runnable locally and emits curated gene-set coverage, full-cell score summaries, donor module scores, primary/strict unadjusted and covariate-adjusted age associations, and module-to-ORA correlations.
- The regulatory-driver path is runnable locally and emits a feasibility note, driver target coverage, pseudobulk and donor-lineage activity scores, age associations, ORA correlations, ranked driver map, and figure.
- The niche-signaling path is runnable locally and emits curated ligand-receptor coverage, donor sender-receiver scores, age associations, ORA correlations, a ranked niche-driver priority table, feasibility note, and figure.

## Exploratory Claims

- Module augmentation slightly improves some ORA model families, but intervals overlap; modules alone are weak in matched feature-family ablation and should be described as supportive biological context rather than the main predictive signal.
- GSE184117 presbyosmia samples show descriptive sample-level shifts in curated modules, marker-only coarse composition panels, marker-reference mapped donor features, and Gateway scANVI/scArches-mapped donor features. The evidence ledger gates these as small-n because the cohort is only 3 versus 3 and no separate author per-cell manual annotation table is available; scANVI mapping confidence is high, but feature-direction concordance is mixed.
- AD/PD donors project to negative ORAA across model families, but this is hypothesis-generating only because each disease group has 5 donors and all disease donors share FLEX v2/device context.
- Genome-wide disease DE has method-sensitive significant rows; limma-voom is more conservative than edgeR and matched analyses reduce sentinel artifacts, but disease biology still needs audited cross-method support.
- A chunked reduced all-cell Gateway scVI run successfully trained on 4,028,275 cells with a 3,000-gene HVG/marker feature set and wrote finite 10-dimensional `X_scvi` coordinates for all cells. Validation passed fine/coarse label-purity and provides the manuscript-primary latent substrate. Cross-run comparison against the 250k seed models and 100k lineage model supports basal and mature-OSN marker continuity, but keeps immature-OSN, progenitor, immune, and sustentacular latent-mechanism wording guarded because top marker labels differ or remain limited across embeddings.
- Full 4M Milo-style neighborhood analyses identify age-associated latent neighborhoods after donor-level adjustment, especially negative age associations in immature OSN/INP/regenerative neighborhoods. Matched FLEX v2/device sensitivity sharply reduces the number of significant neighborhoods but preserves a negative Early iOSN/iOSN signal, while secretory-only neighborhoods do not survive matched correction. ORA-theme annotation separates matched regenerative-neuronal support from all-donor hypothesis-map themes and matched immune-support themes. Curated program scoring shows the matched significant Early iOSN neighborhood is enriched for immature-neuron genes and depleted for HBC programs. Age-bin robustness and exact-neighborhood edgeR count-model parity support the direction of the Python Milo-style result, including the matched Early iOSN neighborhood. Official MiloR subset sensitivity confirms broad age-associated lineage-neighborhood structure but does not independently make matched Early iOSN the dominant signal, so the neighborhood result remains a conservative secondary mechanistic layer rather than a primary discovery claim.

## Deferred Claims

- Larger donor-level external validation of ORA-associated cell-state composition features in independent olfactory aging or NDD datasets.
- Pseudotime, lineage-density, and cNMF claims until method-specific workflows and sensitivity checks pass. Milo-style neighborhood results are available from the full 4M run, but official MiloR sensitivity narrows the interpretation; do not promote exact Early iOSN depletion as an official-MiloR discovery.
- Causal pathway-driver claims from regeneration modules alone; the module layer nominates mechanisms but does not prove signaling activity, lineage flux, or perturbation response.
- Treating the first-pass regulatory-driver layer as SCENIC/decoupler regulon inference; it is a curated target-program fallback until external-prior regulon methods are installed, run, and sensitivity-checked.
- Treating curated ligand-receptor pseudobulk scores as full LIANA/NicheNet/CellChat/CellPhoneDB inference, physical proximity, or spatially resolved signaling; this layer is a transparent fallback until full LR/background testing or spatial/perturbation validation is completed.
- Disease-specific biological interpretation from genome-wide DE without matched/sensitivity audits and edgeR/limma parity.

## Prohibited Framing For Current State

- Do not call ORA an absolute biological-age clock, exact biological-age estimator, or universal aging clock.
- Do not claim context-invariant performance across collection methods, sex groups, ethnicity strata, or sites.
- Do not claim full olfactory specificity until airway/lung/nasal, skin/gut, and immune comparator age-effect estimates are actually computed.
- Do not claim AD/PD diagnostic utility, prediction, classification, biomarker status, or disease specificity.
- Do not claim that cross-sectional donor composition measures lineage flux, cell fate conversion, or regeneration rate.
- Do not interpret UMAP-only structure as trajectory, pseudotime, or lineage evidence.
- Do not call GSE184117 an independent author-label replication of ORA; the authors confirmed no separate per-cell manual annotation table is available. Keep it framed as public-matrix reanalysis/contextual support unless future independent labels or larger donor cohorts are recovered.
- Do not promote exact Early iOSN depletion as a standalone official-MiloR discovery; keep it framed as a narrow exact-neighborhood subclaim with separate parity evidence.
- Do not call niche ligand-receptor scores causal, spatial, or activated signaling mechanisms; keep them framed as prioritized hypotheses.

## Final Claim-Audit Checklist

- [ ] Abstract states healthy-donor ORA as modest, reproducible, and interpretable; it does not use clock, diagnostic, or causal regeneration wording.
- [ ] Results separate the primary healthy-aging ORA claim from external validation, AD/PD projection, DE, and latent-neighborhood extensions.
- [ ] Discussion repeats that ORA is a relative tissue-state axis and lists the main limitations: under-dispersion, cross-sectional design, small external validation, small/confounded disease cohorts, and guarded latent-neighborhood interpretation.
- [ ] Main and extended figure captions avoid diagnostic, trajectory, fate, lineage-flux, and independent-replication wording.
- [ ] Supplementary tables preserve claim-gate, validation-strength, backend/fallback, and provenance columns where applicable.
- [ ] Any new text that mentions AD/PD, GSE184117, Milo/MiloR, scVI, UMAP, pseudotime, or lineage remodeling is checked against this ledger before submission.
