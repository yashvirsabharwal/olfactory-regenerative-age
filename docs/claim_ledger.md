# ORA Claim Ledger

Updated: 2026-06-22

## Supported Primary Claim

- Healthy donor olfactory epithelial composition contains a reproducible, modest age-associated regenerative-state signal above shuffled-age null models.
- Stable ORA features can be organized into supporting/secretory epithelial, immune/inflammatory, stress/senescence, neuronal-lineage, and regenerative/progenitor themes, with each row treated as associational rather than causal.

## Supported Engineering Claim

- The repository can generate donor-level features, ORA model benchmarks, module summaries, NDD projections, targeted pseudobulk DE, genome-wide pseudobulk export, edgeR and limma-voom summaries, and an MVP from ignored local artifacts.

## Exploratory Claims

- Module augmentation slightly improves some ORA model families, but intervals overlap.
- GSE184117 presbyosmia samples show descriptive sample-level shifts in curated modules, marker-only coarse composition panels, marker-reference mapped donor features, and Gateway scANVI/scArches-mapped donor features. The evidence ledger gates these as small-n because the cohort is only 3 versus 3; scANVI mapping confidence is high, but feature-direction concordance is mixed.
- AD/PD donors project to negative ORAA across model families, but this is hypothesis-generating only because each disease group has 5 donors and all disease donors share FLEX v2/device context.
- Genome-wide disease DE has method-sensitive significant rows; limma-voom is more conservative than edgeR and matched analyses reduce sentinel artifacts, but disease biology still needs audited cross-method support.
- A chunked reduced all-cell Gateway scVI run successfully trained on 4,028,275 cells with a 3,000-gene HVG/marker feature set and wrote finite 10-dimensional `X_scvi` coordinates for all cells. Validation passed fine/coarse label-purity and provides the manuscript-primary latent substrate. Cross-run comparison against the 250k seed models and 100k lineage model supports basal and mature-OSN marker continuity, but keeps immature-OSN, progenitor, immune, and sustentacular latent-mechanism wording guarded because top marker labels differ or remain limited across embeddings.
- Full 4M Milo-style neighborhood analyses identify age-associated latent neighborhoods after donor-level adjustment, especially negative age associations in immature OSN/INP/regenerative neighborhoods. Matched FLEX v2/device sensitivity sharply reduces the number of significant neighborhoods but preserves a negative Early iOSN/iOSN signal, while secretory-only neighborhoods do not survive matched correction. ORA-theme annotation separates matched regenerative-neuronal support from all-donor hypothesis-map themes and matched immune-support themes. Curated program scoring shows the matched significant Early iOSN neighborhood is enriched for immature-neuron genes and depleted for HBC programs. Age-bin robustness and exact-neighborhood edgeR count-model parity support the direction of the Python Milo-style result, including the matched Early iOSN neighborhood. Official MiloR subset sensitivity confirms broad age-associated lineage-neighborhood structure but does not independently make matched Early iOSN the dominant signal, so the neighborhood result remains a conservative secondary mechanistic layer rather than a primary discovery claim.

## Deferred Claims

- Larger donor-level external validation of ORA-associated cell-state composition features in independent olfactory aging or NDD datasets.
- Pseudotime, lineage-density, and cNMF claims until method-specific workflows and sensitivity checks pass. Milo-style neighborhood results are available from the full 4M run, but official MiloR sensitivity narrows the interpretation; do not promote exact Early iOSN depletion as an official-MiloR discovery.
- Disease-specific biological interpretation from genome-wide DE without matched/sensitivity audits and edgeR/limma parity.

## Prohibited Framing For Current State

- Do not call ORA an absolute biological-age clock.
- Do not claim AD/PD diagnostic utility.
- Do not claim measured lineage flux from cross-sectional composition.
- Do not interpret UMAP-only structure as trajectory evidence.
