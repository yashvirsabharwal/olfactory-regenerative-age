# scVI Embedding Comparison Claim Gate

Updated: 2026-06-22

## Verdict

The primary full 4M reduced scVI model represents 4,028,275 cells, has 10 latent dimensions, fine-label purity 0.877, and coarse-label purity 0.975.

The full 4M reduced model remains the manuscript-primary latent substrate. The 250k seed models and the 100k lineage model should be used as sensitivity anchors, not as competing primary analyses.

Supported marker panels across the current comparison: basal, mature_osn.
Marker-specific caveats remain for: immature_osn, immune, progenitor, sustentacular.

## Manuscript Rule

- Use full 4M reduced scVI for the main latent/neighborhood methods.
- Use 250k seed and lineage-focused runs to qualify marker-continuity and lineage-specific interpretation.
- Keep progenitor and immune latent-mechanism claims guarded unless follow-up marker/program or compartment-specific analyses support them.
- Describe the Early iOSN neighborhood as an exact-neighborhood finding supported by edgeR, age-bin directionality, and curated program scoring, not as a fully replicated official-MiloR or all-embedding claim.
