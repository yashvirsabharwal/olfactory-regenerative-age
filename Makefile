PYTHON ?= python3

R_ENV ?= .mamba/ora-r
MICROMAMBA ?= $(HOME)/.local/bin/micromamba
RSCRIPT ?= $(MICROMAMBA) run -p $(R_ENV) Rscript

.PHONY: setup test download-gateway download-info toy-data smoke-toy inspect cohort aggregate features features-augmented age-associations model-ora model-ora-repeated model-ora-augmented project-ndd project-ndd-uncertainty report modules published-gene-modules external-validation trajectory pseudobulk pseudobulk-genomewide pseudobulk-genomewide-qc pseudobulk-genomewide-edger pseudobulk-genomewide-de-summary pseudobulk-covariate-de ora-sensitivity milo cnmf clean

setup:
	$(PYTHON) -m pip install -e ".[dev]"

test:
	$(PYTHON) -m unittest discover -s tests -p "test_*.py"

download-gateway:
	$(PYTHON) scripts/download_cellxgene.py --config configs/gateway.yaml

download-info:
	$(PYTHON) scripts/download_cellxgene.py --config configs/gateway.yaml --dry-run

toy-data:
	$(PYTHON) scripts/create_toy_gateway_h5ad.py --out data/raw/toy_gateway.h5ad

smoke-toy:
	$(PYTHON) scripts/run_toy_smoke.py

inspect:
	$(PYTHON) scripts/inspect_h5ad.py --config configs/gateway.yaml

cohort:
	$(PYTHON) scripts/build_sample_manifest.py --config configs/gateway.yaml

aggregate:
	$(PYTHON) scripts/aggregate_cell_counts.py --config configs/gateway.yaml

features:
	$(PYTHON) scripts/build_feature_matrix.py --config configs/gateway.yaml

features-augmented:
	$(PYTHON) scripts/build_feature_matrix.py --config configs/gateway.yaml --include-modules

age-associations:
	$(PYTHON) scripts/run_age_associations.py --features data/processed/donor_cell_state_features.tsv --manifest data/processed/cohort_manifest.tsv --config configs/gateway.yaml

model-ora:
	$(PYTHON) scripts/run_age_models.py --features data/processed/ora_feature_matrix.tsv --manifest data/processed/cohort_manifest.tsv --config configs/models.yaml

model-ora-repeated:
	$(PYTHON) scripts/run_age_models_repeated.py --features data/processed/ora_feature_matrix.tsv --manifest data/processed/cohort_manifest.tsv --model-config configs/models.yaml --gateway-config configs/gateway.yaml

model-ora-augmented:
	$(PYTHON) scripts/run_age_models.py --features data/processed/ora_augmented_feature_matrix.tsv --manifest data/processed/cohort_manifest.tsv --config configs/models.yaml --performance-out results/tables/ora_augmented_model_performance.tsv --scores-out results/tables/augmented_donor_ora_scores.tsv --importance-out results/tables/ora_augmented_feature_importance.tsv

project-ndd:
	$(PYTHON) scripts/project_ndd_ora.py --gateway-config configs/gateway.yaml --model-config configs/models.yaml

project-ndd-uncertainty:
	$(PYTHON) scripts/summarize_ndd_projection_uncertainty.py --gateway-config configs/gateway.yaml

report:
	$(PYTHON) scripts/generate_mvp_report.py --config configs/gateway.yaml

trajectory:
	$(PYTHON) scripts/run_trajectory.py --config configs/gateway.yaml

modules:
	$(PYTHON) scripts/score_gene_sets.py --config configs/gateway.yaml --gene-sets configs/gene_sets.yaml

published-gene-modules:
	$(PYTHON) scripts/score_gene_sets.py --config configs/gateway.yaml --gene-sets configs/external_datasets.yaml --summary-out results/tables/published_module_score_summary.tsv --coverage-out results/tables/published_module_gene_list_coverage.tsv --donor-features-out data/processed/published_donor_module_features.tsv

external-validation:
	$(PYTHON) scripts/summarize_external_validation.py --external-config configs/external_datasets.yaml --gateway-config configs/gateway.yaml

pseudobulk:
	$(PYTHON) scripts/aggregate_pseudobulk.py --config configs/gateway.yaml --gene-sets configs/gene_sets.yaml

pseudobulk-genomewide:
	$(PYTHON) scripts/export_genomewide_pseudobulk.py --config configs/gateway.yaml

pseudobulk-genomewide-qc:
	$(PYTHON) scripts/summarize_genomewide_pseudobulk.py --config configs/gateway.yaml

pseudobulk-genomewide-edger:
	$(RSCRIPT) scripts/run_genomewide_edger.R --counts data/processed/pseudobulk_genomewide_counts.tsv.gz --metadata data/processed/pseudobulk_genomewide_metadata.tsv --manifest data/processed/cohort_manifest.tsv --out results/tables/pseudobulk_genomewide_edger.tsv.gz --summary-out results/tables/pseudobulk_genomewide_edger_summary.tsv

pseudobulk-genomewide-de-summary:
	$(PYTHON) scripts/summarize_genomewide_de.py --config configs/gateway.yaml

pseudobulk-covariate-de:
	$(PYTHON) scripts/run_pseudobulk_covariate_de.py --config configs/gateway.yaml

ora-sensitivity:
	$(PYTHON) scripts/run_ora_sensitivity.py --gateway-config configs/gateway.yaml --model-config configs/models.yaml

milo:
	Rscript scripts/run_milo.R configs/gateway.yaml

cnmf:
	$(PYTHON) scripts/run_cnmf.py --config configs/gateway.yaml

clean:
	rm -rf .pytest_cache .ruff_cache build dist *.egg-info
