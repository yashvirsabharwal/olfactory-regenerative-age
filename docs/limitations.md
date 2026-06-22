# ORA Limitations Draft

Updated: 2026-06-22

- The analysis is cross-sectional and cannot measure lineage flux directly.
- ORA predictions are under-dispersed, so ORA should be interpreted as a relative regenerative-state axis rather than an absolute age estimator.
- AD and PD projection cohorts have 5 donors each and are fully tied to FLEX v2/device collection context.
- Genome-wide DE is sensitive to sex, chemistry, collection method, donor-balance differences, and statistical engine; top hits require audit, matched/sensitivity checks, and edgeR/limma agreement before biological emphasis.
- DESeq2 is installed but deferred for the full fine-cell-state disease-DE grid; use it only for larger validation cohorts or targeted shortlisted hits with manually reviewed diagnostics.
- External validation is small-n; `GSE184117` now supports raw sample-level module scoring, marker-only composition, marker-reference mapped donor features, Gateway scANVI/scArches donor features, and mapped-feature concordance, but the cohort is still only 3 versus 3 and scANVI feature concordance is mixed.
- The CELLxGENE export currently exposes UMAP but not author scVI/scANVI latent coordinates. A chunked reduced all-cell scVI run now provides local `X_scvi` coordinates for all Gateway cells, and comparison against 250k seed and 100k lineage models keeps the full 4M model as manuscript-primary. Trajectory claims remain deferred, and immature-OSN/progenitor/immune/sustentacular latent-mechanism wording remains guarded where marker continuity differs across embeddings.
- Module scores use average log1p expression over curated genes, not exact rank-based UCell scores.
