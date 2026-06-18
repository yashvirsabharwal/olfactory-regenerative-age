# ORA Claim Ledger

Updated: 2026-06-18

## Supported Primary Claim

- Healthy donor olfactory epithelial composition contains a reproducible, modest age-associated regenerative-state signal above shuffled-age null models.
- Stable ORA features can be organized into supporting/secretory epithelial, immune/inflammatory, stress/senescence, neuronal-lineage, and regenerative/progenitor themes, with each row treated as associational rather than causal.

## Supported Engineering Claim

- The repository can generate donor-level features, ORA model benchmarks, module summaries, NDD projections, targeted pseudobulk DE, genome-wide pseudobulk export, edgeR and limma-voom summaries, and an MVP report from ignored local artifacts.

## Exploratory Claims

- Module augmentation slightly improves some ORA model families, but intervals overlap.
- GSE184117 presbyosmia samples show descriptive sample-level shifts in curated modules, marker-only coarse composition panels, marker-reference mapped donor features, and Gateway scANVI/scArches-mapped donor features. The evidence ledger gates these as small-n because the cohort is only 3 versus 3; scANVI mapping confidence is high, but feature-direction concordance is mixed.
- AD/PD donors project to negative ORAA across model families, but this is hypothesis-generating only because each disease group has 5 donors and all disease donors share FLEX v2/device context.
- Genome-wide disease DE has method-sensitive significant rows; limma-voom is more conservative than edgeR and matched analyses reduce sentinel artifacts, but disease biology still needs audited cross-method support.
- A chunked reduced all-cell Gateway scVI run successfully trained on 4,028,275 cells with a 3,000-gene HVG/marker feature set and wrote finite 10-dimensional `X_scvi` coordinates for all cells. Validation passed fine/coarse label-purity, FLEX/device/condition/sex mixing, and basal/OSN/sustentacular marker-continuity checks; progenitor and immune continuity remain limited in the full model and should be cross-checked against the 250k seed and lineage-focused models before mechanistic claims.
- Full 4M Milo-style neighborhood analyses identify age-associated latent neighborhoods after donor-level adjustment, especially negative age associations in immature OSN/INP/regenerative neighborhoods. Matched FLEX v2/device sensitivity sharply reduces the number of significant neighborhoods but preserves a negative Early iOSN/iOSN signal, while secretory-only neighborhoods do not survive matched correction. This supports a conservative exploratory mechanistic layer until neighborhood marker annotation and implementation-parity checks are complete.

## Deferred Claims

- Larger donor-level external validation of ORA-associated cell-state composition features in independent olfactory aging or NDD datasets.
- Pseudotime, lineage-density, and cNMF claims until method-specific workflows and sensitivity checks pass. Milo-style neighborhood results are available from the full 4M run but remain exploratory until marker/program annotation and implementation-parity checks pass.
- Disease-specific biological interpretation from genome-wide DE without matched/sensitivity audits and edgeR/limma parity.

## Prohibited Framing For Current State

- Do not call ORA an absolute biological-age clock.
- Do not claim AD/PD diagnostic utility.
- Do not claim measured lineage flux from cross-sectional composition.
- Do not interpret UMAP-only structure as trajectory evidence.
