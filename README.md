# Olfactory Regenerative Age

This repository is a publication-oriented donor-level analysis of olfactory regenerative aging in the Gateway human olfactory epithelium atlas.

Primary source: **Gateway: patient olfactory neurons for large-scale discovery in neurodegenerative disease**, DOI `10.64898/2026.06.10.731272`.

## Current Scope

The implementation focuses on memory-safe metadata inspection, healthy-donor cohort definition, cell-state composition features, age associations, composition-only and module-augmented ORA modeling, frozen healthy-trained AD/PD projection, reporting, chunked average-expression module scoring, targeted curated-gene pseudobulk DE, donor-level covariate-adjusted pseudobulk DE, genome-wide edgeR and limma-voom DE summaries, GSE184117 raw 10x mapping, Gateway scANVI/scArches reference mapping, scaled and all-cell reduced scVI latent validation, full-scale Milo-style neighborhood analyses with edgeR and official-MiloR subset sensitivity, provenance, and manuscript-readiness tables. Deferred mechanistic extensions such as pseudotime, CellRank, and cNMF are tracked in `docs/active_work_tracker.md`.

Manuscript framing is tracked in `docs/manuscript_framework.md`; journal-readiness gates are tracked in `docs/journal_acceptance_tracker.md`. External-validation exhaustion and reproducibility packaging are recorded in `docs/external_validation_final_search.md`, `docs/gse184117_label_request.md`, `docs/large_artifact_manifest.md`, and `docs/manuscript_rerun_profile.md`.

## Quickstart

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[dev]"
make test
```

The Makefile automatically prefers `.venv/bin/python` when it exists. If your shell is unusual, the explicit form is:

```bash
PYTHON=.venv/bin/python make test
```

Place the Gateway H5AD at `data/raw/gateway.h5ad`, or update `configs/gateway.yaml`.

```bash
make download-gateway
make inspect
make cohort
make aggregate
make features
make age-associations
make model-ora
make modules
make features-augmented
make model-ora-augmented
make project-ndd
make pseudobulk
make pseudobulk-covariate-de
make pseudobulk-genomewide-de-summary
make pseudobulk-genomewide-de-audit
make pseudobulk-genomewide-edger-matched
make pseudobulk-genomewide-de-summary-matched
make pseudobulk-genomewide-de-audit-matched
make pseudobulk-genomewide-limma
make pseudobulk-genomewide-limma-de-summary
make pseudobulk-genomewide-limma-de-audit
make pseudobulk-genomewide-limma-matched
make pseudobulk-genomewide-limma-de-summary-matched
make pseudobulk-genomewide-limma-de-audit-matched
make external-validation
make external-gse184117-modules
make external-gse184117-markers
make external-marker-age-concordance
make external-evidence
make feature-interpretation
make latent-space-audit
make latent-space-recompute-plan
make scvi-scaled-250k
make scvi-scaled-validation
make scvi-reduced-4m
make scvi-full-4m-reduced
make scvi-lineage-basal-neural
make scvi-lineage-validation
make model-card
make manuscript-figures
make publication-tables
make output-provenance
make report
```

If you do not have the real Gateway file yet, run the smoke workflow against a temporary toy H5AD:

```bash
make smoke-toy
```

`make download-gateway` uses the resolved CELLxGENE H5AD URL in `configs/gateway.yaml`. The file is large: 27,651,770,343 bytes, about 25.8 GiB. To check the resolved source without downloading:

```bash
make download-info
```

To download and inventory the first external aging/presbyosmia validation archive:

```bash
make download-gse184117
make inspect-gse184117
make external-gse184117-modules
make external-gse184117-markers
make external-gse184117-mapped
make external-marker-age-concordance
make external-mapped-feature-concordance
make external-scanvi-reference-map
make external-scanvi-feature-concordance
make external-validation
make external-evidence
```

## Outputs

Generated data and analysis outputs are intentionally ignored by Git:

- `data/processed/cohort_manifest.tsv`
- `data/processed/donor_cell_state_counts.tsv`
- `data/processed/donor_cell_state_features.tsv`
- `data/processed/donor_module_features.tsv`
- `data/processed/ora_feature_matrix.tsv`
- `data/processed/ora_augmented_feature_matrix.tsv`
- `results/tables/age_cell_state_associations.tsv`
- `results/tables/ora_model_performance.tsv`
- `results/tables/donor_ora_scores.tsv`
- `results/tables/ora_augmented_model_performance.tsv`
- `results/tables/augmented_donor_ora_scores.tsv`
- `results/tables/ora_augmented_feature_importance.tsv`
- `results/tables/ora_calibration.tsv`
- `results/tables/ora_calibrated_scores.tsv`
- `results/tables/ora_augmented_candidate_repeated_cv_summary.tsv`
- `results/tables/ndd_ora_projection.tsv`
- `results/tables/ndd_ora_projection_summary.tsv`
- `results/tables/module_score_summary.tsv`
- `results/tables/module_gene_coverage.tsv`
- `data/processed/pseudobulk_counts.tsv.gz`
- `data/processed/pseudobulk_metadata.tsv`
- `results/tables/pseudobulk_gene_coverage.tsv`
- `results/tables/pseudobulk_de.tsv`
- `results/tables/pseudobulk_covariate_de.tsv`
- `results/tables/pseudobulk_genomewide_de_audit.tsv`
- `results/tables/pseudobulk_genomewide_donor_balance.tsv`
- `results/tables/pseudobulk_genomewide_matched_feasibility.tsv`
- `results/tables/pseudobulk_genomewide_edger_matched_flex_v2_device.tsv.gz`
- `results/tables/pseudobulk_genomewide_de_summary_matched_flex_v2_device.tsv`
- `results/tables/pseudobulk_genomewide_de_audit_matched_flex_v2_device.tsv`
- `results/tables/pseudobulk_genomewide_limma_voom.tsv.gz`
- `results/tables/pseudobulk_genomewide_limma_voom_de_summary.tsv`
- `results/tables/pseudobulk_genomewide_limma_voom_de_audit.tsv`
- `results/tables/pseudobulk_genomewide_limma_voom_matched_flex_v2_device.tsv.gz`
- `results/tables/pseudobulk_genomewide_limma_voom_de_summary_matched_flex_v2_device.tsv`
- `results/tables/pseudobulk_genomewide_limma_voom_de_audit_matched_flex_v2_device.tsv`
- `results/tables/ora_model_card.tsv`
- `results/tables/ndd_label_permutation.tsv`
- `results/tables/external_raw_inventory.tsv`
- `results/tables/external_sample_metadata.tsv`
- `results/tables/external_10x_sample_qc.tsv`
- `results/tables/external_10x_module_scores.tsv`
- `results/tables/external_10x_module_contrasts.tsv`
- `results/tables/external_10x_marker_coverage.tsv`
- `results/tables/external_10x_marker_composition.tsv`
- `results/tables/external_10x_marker_contrasts.tsv`
- `results/tables/external_marker_age_concordance.tsv`
- `data/processed/gse184117_marker_mapped.h5ad`
- `results/tables/external_10x_mapping_qc.tsv`
- `data/processed/gse184117_mapped_donor_features.tsv`
- `results/tables/external_mapped_feature_concordance.tsv`
- `data/processed/gse184117_scanvi_mapped.h5ad`
- `results/tables/external_scanvi_mapping_qc.tsv`
- `data/processed/gse184117_scanvi_donor_features.tsv`
- `results/tables/external_scanvi_feature_concordance.tsv`
- `results/tables/external_validation_evidence.tsv`
- `results/tables/ora_feature_interpretation.tsv`
- `results/tables/latent_space_local_audit.tsv`
- `results/tables/latent_space_portal_assets.tsv`
- `results/tables/latent_space_readiness.tsv`
- `results/tables/latent_recompute_feasibility.tsv`
- `results/tables/scvi_pilot_validation.tsv`
- `data/processed/gateway_scvi_stratified_250k.h5ad`
- `results/tables/scvi_scaled_250k_validation.tsv`
- `data/processed/gateway_hvg3003_4m.h5ad`
- `data/processed/gateway_scvi_full_4m_reduced.h5ad`
- `results/tables/scvi_full_4m_reduced_validation.tsv`
- `results/tables/scvi_scaled_250k_seed23_validation.tsv`
- `results/tables/scvi_embedding_claim_gates.tsv`
- `results/tables/scvi_embedding_marker_concordance.tsv`
- `data/processed/gateway_scvi_lineage_basal_neural_100k.h5ad`
- `results/tables/scvi_lineage_basal_neural_validation.tsv`
- `docs/latent_space_recovery_plan.md`
- `docs/latent_recompute_workflow.md`
- `results/reports/command_manifest.tsv`
- `results/reports/output_provenance.tsv`
- `results/reports/mvp_report.md`
- `results/figures/mvp_*.png`
- `results/figures/manuscript_figure*.png`
- `results/figures/manuscript_figure*.pdf`
- `results/figures/extended_data_figure*.png`
- `results/figures/extended_data_figure*.pdf`
- `results/tables/manuscript_table_*.tsv`
- `docs/publication_tables.md`
- `docs/manuscript_framework.md`
- `manuscript/main.tex`
- `manuscript/references.bib`

## Guardrails

- Train ORA only on healthy donors with valid age.
- Keep AD/PD donors for frozen-model projection only.
- Split and cross-validate by donor, never by cell.
- Treat trajectory output as a cross-sectional lineage-density proxy, not measured lineage flux.
- Use chemistry, collection method, site, and yield as covariates or sensitivity variables, not primary biological ORA features.
