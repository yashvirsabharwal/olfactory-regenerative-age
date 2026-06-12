PYTHON ?= python3

.PHONY: setup test download-gateway download-info toy-data smoke-toy inspect cohort aggregate features features-augmented age-associations model-ora model-ora-augmented report modules trajectory pseudobulk milo cnmf clean

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

model-ora-augmented:
	$(PYTHON) scripts/run_age_models.py --features data/processed/ora_augmented_feature_matrix.tsv --manifest data/processed/cohort_manifest.tsv --config configs/models.yaml --performance-out results/tables/ora_augmented_model_performance.tsv --scores-out results/tables/augmented_donor_ora_scores.tsv --importance-out results/tables/ora_augmented_feature_importance.tsv

report:
	$(PYTHON) scripts/generate_mvp_report.py --config configs/gateway.yaml

trajectory:
	$(PYTHON) scripts/run_trajectory.py --config configs/gateway.yaml

modules:
	$(PYTHON) scripts/score_gene_sets.py --config configs/gateway.yaml --gene-sets configs/gene_sets.yaml

pseudobulk:
	$(PYTHON) scripts/aggregate_pseudobulk.py --config configs/gateway.yaml --gene-sets configs/gene_sets.yaml

milo:
	Rscript scripts/run_milo.R configs/gateway.yaml

cnmf:
	$(PYTHON) scripts/run_cnmf.py --config configs/gateway.yaml

clean:
	rm -rf .pytest_cache .ruff_cache build dist *.egg-info
