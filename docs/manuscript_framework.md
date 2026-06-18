# ORA Manuscript Framework

Updated: 2026-06-18

## Working Title

An interpretable olfactory regenerative aging axis from human single-cell epithelial composition

## Target Article Type

Primary target: computational single-cell atlas reanalysis/resource manuscript with a healthy-aging biology centerpiece.

The paper should be framed as a reproducible and interpretable donor-level analysis of human olfactory epithelial aging. It should not be framed as a clinical age clock, a diagnostic neurodegeneration classifier, or a trajectory/neighborhood-discovery paper in the current state.

## Main Text Structure

1. Abstract
   - One-sentence problem: the aging human olfactory epithelium is clinically important but donor-level regenerative-state readouts are underdeveloped.
   - One-sentence approach: reanalyze the Gateway 4M-cell atlas into donor-level composition, lineage-ratio, and curated module features.
   - One-sentence result: healthy-donor ORA models capture a modest but reproducible aging signal above shuffled-age nulls, with interpretable features.
   - One-sentence guardrail: NDD projection, external validation, genome-wide DE, and latent-space analyses are explicitly claim-gated.

2. Introduction
   - Human olfactory epithelium as a renewing neuroepithelium and aging-relevant sensory barrier.
   - Why donor-level single-cell composition is a useful complement to bulk clocks and molecular age predictors.
   - Gap: existing olfactory datasets are rich, but regenerative aging has not been formalized as an interpretable donor-level axis.
   - Contribution: ORA as a reproducible, audited framework for healthy olfactory epithelial aging.

3. Results
   - Result 1: Gateway cohort and analysis design.
     - 202 donors, 4,028,275 cells, 187 healthy age-usable training donors.
     - Healthy-only model training; AD/PD held out for frozen projection.
   - Result 2: Donor-level cell-state composition captures age-associated regenerative structure.
     - Age-associated proportions, CLR features, and lineage ratios.
     - Biological themes: supporting/secretory epithelium, immune/inflammatory state, stress/senescence, neuronal maturation, basal/progenitor state.
   - Result 3: ORA models learn a modest but reproducible aging axis.
     - Repeated CV, shuffled-age null, nested tuning, calibration diagnostics.
     - Emphasize relative regenerative-state axis, not absolute chronological-age clock.
   - Result 4: Curated modules add biological annotation and only modest predictive gain.
     - Module coverage and module-augmented benchmarks.
     - Feature interpretation table links stable predictors to tissue biology.
   - Result 5: External olfactory aging evidence is supportive but not definitive.
     - GSE184117 raw 10x sample-level module, marker-only composition, marker-reference mapped donor features, scANVI/scArches donor features, and mapped-feature concordance.
     - Marker-age, marker-reference mapped-feature, and scANVI-mapped feature concordance against Gateway age features.
     - Explicit small-n/non-scANVI-transfer claim gate.
   - Result 6: AD/PD ORA projection is directionally stable but exploratory.
     - Frozen healthy-trained projections only.
     - Matched FLEX v2/device sensitivity and label permutation guardrails.
   - Result 7: Genome-wide pseudobulk disease DE is audited hypothesis generation.
     - edgeR and limma-voom parity.
     - Donor balance, sex-linked sentinel, matched subset, and method sensitivity.
   - Result 8: Full-scale latent neighborhoods reveal an exploratory regenerative aging layer.
     - CELLxGENE export has only `X_umap`.
     - A chunked reduced all-cell scVI run produced finite 10-dimensional `X_scvi` for all 4,028,275 Gateway cells, with 250k seed and lineage-focused models as sensitivity anchors.
     - Full 4M Milo-style neighborhood analyses identify age-associated neighborhoods, especially negative age associations in immature OSN/INP/regenerative neighborhoods.
     - Matched technical sensitivity, marker/program annotation, and implementation parity remain required before promoting the neighborhood result as a primary claim.

4. Discussion
   - Interpretation: ORA detects a tissue-specific regenerative aging axis in a renewing human neuroepithelium.
   - Novelty: interpretable donor-level single-cell composition rather than black-box universal age prediction.
   - Why modest accuracy is biologically meaningful: under-dispersed predictions and calibration caveats reflect tissue-state signal rather than exact age reconstruction.
   - External validation and disease extension remain the key upgrades.
   - Practical resource value: command manifest, model card, claim ledger, audited DE, and validation-strength reporting.

5. Methods
   - Dataset acquisition and metadata harmonization.
   - Cohort definitions and donor-level splitting.
   - Donor-level composition, CLR, lineage-ratio, and module features.
   - ORA model families, repeated CV, nested tuning, null testing, calibration, residual diagnostics.
   - External GSE184117 raw-adapter, marker-reference mapped donor features, scANVI/scArches query mapping, and mapped-feature concordance.
   - NDD frozen projection and label-permutation diagnostics.
   - Pseudobulk construction, edgeR QL, limma-voom parity, matched sensitivity, and audit tables.
   - Latent-space audit, scaled scVI validation, lineage scVI validation, and deferred trajectory gates.
   - Software, command provenance, and artifact policy.

6. Limitations
   - Cross-sectional analysis cannot measure lineage flux.
   - ORA is under-dispersed and should be interpreted as a relative axis.
   - External validation is currently small-n; scANVI/scArches transfer is available, but feature-direction concordance is mixed.
   - AD/PD cohorts are very small and technically confounded.
   - Genome-wide DE is hypothesis-generating until replicated and audited.
   - Latent-space analyses are scaled-QC stage but not yet seed-stable mechanism claims.

7. Data And Code Availability
   - Gateway source and CELLxGENE identifiers.
   - GEO accessions for GSE184117 and GSE151973.
   - Repository branch, command manifest, output provenance, ignored large artifacts.

## Main Figures

1. Figure 1: cohort, workflow, and claim-gated analysis design.
2. Figure 2: age-associated donor-level composition and interpretable tissue-state themes.
3. Figure 3: ORA model performance, null testing, and calibration.
4. Figure 4: stable feature interpretation and module-augmented biology.
5. Figure 5: external validation and NDD projection guardrails.
6. Figure 6: genome-wide DE audit and latent-space readiness.

## Extended Data

- Extended Data 1: full model-card table.
- Extended Data 2: residual diagnostics by age bin, sex, chemistry, collection method, race/ethnicity, site, and yield.
- Extended Data 3: full external validation evidence ledger, marker-age concordance, marker-reference mapped-feature concordance, and scANVI-mapped feature concordance.
- Extended Data 4: NDD donor appendix and matched FLEX v2/device projection sensitivity.
- Extended Data 5: edgeR and limma-voom DE top-hit and audit tables.
- Extended Data 6: all-cell, scaled-seed, and lineage scVI validation tables plus full 4M Milo-style neighborhood maps and claim gates.

## Current Publication Standing

The project is on track for a conservative manuscript after figure polishing and LaTeX drafting. The most defensible paper is already substantial: a large human olfactory atlas reanalysis with interpretable healthy-aging signal, robust model benchmarking, explicit null testing, external sanity checks, and audited disease extensions. The manuscript should lead with the healthy aging axis and use external validation, NDD projection, DE, and scVI as guarded support rather than central claims.
