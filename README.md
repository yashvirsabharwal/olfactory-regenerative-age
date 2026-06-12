# Olfactory Regenerative Age

This repository is an MVP scaffold for a donor-level analysis of olfactory regenerative aging in the Gateway human olfactory epithelium atlas.

Primary source: **Gateway: patient olfactory neurons for large-scale discovery in neurodegenerative disease**, DOI `10.64898/2026.06.10.731272`.

## MVP Scope

The initial implementation focuses on memory-safe metadata inspection, healthy-donor cohort definition, cell-state composition features, age associations, composition-only ORA modeling, reporting, and chunked average-expression module scoring. Heavier analyses such as exact rank-based UCell scoring, pseudobulk differential expression, pseudotime density, Milo, and cNMF remain exposed as explicit deferred commands.

## Quickstart

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[dev]"
make test
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
- `results/tables/module_score_summary.tsv`
- `results/tables/module_gene_coverage.tsv`
- `results/reports/mvp_report.md`
- `results/figures/mvp_*.png`

## Guardrails

- Train ORA only on healthy donors with valid age.
- Keep AD/PD donors for frozen-model projection only.
- Split and cross-validate by donor, never by cell.
- Treat trajectory output as a cross-sectional lineage-density proxy, not measured lineage flux.
- Use chemistry, collection method, site, and yield as covariates or sensitivity variables, not primary biological ORA features.
