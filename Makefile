PYTHON ?= $(shell if [ -x .venv/bin/python ]; then printf ".venv/bin/python"; else printf "python3"; fi)

R_ENV ?= .mamba/ora-r
MICROMAMBA ?= $(HOME)/.local/bin/micromamba
RSCRIPT ?= $(MICROMAMBA) run -p $(R_ENV) Rscript
MILO_H5AD ?= data/processed/gateway_scvi_stratified_250k.h5ad
MILO_FULL_H5AD ?= data/processed/gateway_scvi_full_4m_reduced.h5ad
MILO_LINEAGE_H5AD ?= data/processed/gateway_scvi_lineage_basal_neural_100k.h5ad
MILO_NEIGHBORHOODS ?= 1000
MILO_FULL_NEIGHBORHOODS ?= 20000
MILO_NEIGHBORS ?= 50
MILO_FULL_NEIGHBORS ?= 100
MILO_MATCHED_DONOR_QUERY ?= chemistry == 'flex_v2' and collection_method == 'device'
MILO_MATCHED_MIN_DONORS ?= 12
MILOR_SUBSET_CELLS ?= 100000
MILOR_MATCHED_SUBSET_CELLS ?= 75000
MILOR_RSCRIPT ?= $(HOME)/bin/micromamba run -p $(HOME)/micromamba-envs/ora-milor Rscript
SCVI_DONOR_H5AD ?= data/processed/gateway_scvi_full_4m_reduced.h5ad
SCVI_DONOR_MODELS ?= hist_gradient_boosting xgboost catboost boosted_ensemble
FOUNDATION_BENCHMARK_H5AD ?= data/raw/gateway.h5ad
FOUNDATION_PYTHON ?= .venv_foundation/bin/python
GENEFORMER_H5AD ?= data/processed/foundation_benchmark_lineage_subset.h5ad
GENEFORMER_MAX_CELLS ?= 24000
GENEFORMER_BATCH_SIZE ?= 64
FOUNDATION_LINEAGE_CELLS ?= 120000
FOUNDATION_EPITHELIAL_CELLS ?= 180000
FOUNDATION_ALLCELL_CELLS ?= 250000
REMOTE_FULL_SCVI_ACTION ?= help

.PHONY: setup test download-gateway download-info download-gse184117 inspect-gse184117 external-gse184117-modules external-gse184117-markers external-gse184117-mapped external-gse184117-status external-scanvi-reference-map external-scanvi-feature-concordance external-mapped-feature-concordance external-marker-age-concordance external-candidate-matrix external-public-data-exhaustion external-evidence inspect cohort aggregate features features-augmented age-associations compositional-age-model negative-controls feature-family-ablation leave-context-out foundation-benchmark-subsets geneformer-donor-features geneformer-age-model foundation-model-benchmark aging-clock-baseline cross-tissue-specificity cross-tissue-age-effects spatial-validation-plan perturbation-validation-plan perturbation-gse309325-organoid regeneration-axis-feature-map regeneration-modules regeneration-module-analysis regulatory-driver-analysis niche-signaling-analysis regeneration-dynamics model-ora model-ora-diagnostics model-ora-repeated model-ora-augmented model-ora-candidate-repeated ora-feature-set-comparison ora-permutation-null ora-nested-tuning ora-stacking feature-interpretation project-ndd project-ndd-feature-sensitivity project-ndd-uncertainty project-ndd-diagnostics project-ndd-label-permutation manuscript manuscript-check manuscript-figures publication-tables modules published-gene-modules external-validation external-feature-harmonization latent-space-audit latent-space-recompute-plan latent-space-plan scvi-pilot scvi-pilot-validation scvi-scaled-250k scvi-scaled-250k-seed23 scvi-scaled-validation scvi-scaled-comparison scvi-embedding-claim-gates scvi-donor-features scvi-donor-age-model scvi-donor-comparison scvi-donor-embedding-baseline scvi-hybrid-features scvi-hybrid-age-model scvi-hybrid-comparison scvi-hybrid-benchmark scvi-scaled-500k scvi-scaled-1m scvi-reduced-4m scvi-full-4m scvi-full-4m-safe scvi-full-4m-reduced scvi-full-validation scvi-lineage-basal-neural scvi-lineage-validation pseudobulk pseudobulk-genomewide pseudobulk-genomewide-qc pseudobulk-genomewide-edger pseudobulk-genomewide-edger-matched pseudobulk-genomewide-limma pseudobulk-genomewide-limma-matched pseudobulk-genomewide-de-summary pseudobulk-genomewide-de-audit pseudobulk-genomewide-de-summary-matched pseudobulk-covariate-de ora-sensitivity ora-sensitivity-rf model-card output-provenance release-manifest environment-report archive-review-package local-light local-medium remote-heavy submission-freeze all-summary remote-full-scvi milo milo-lineage milo-secretory milo-full-4m milo-full-4m-lineage milo-full-4m-secretory milo-full-4m-matched milo-full-4m-lineage-matched milo-full-4m-secretory-matched milo-full-4m-lineage-programs milo-full-4m-lineage-matched-programs milo-full-4m-lineage-age-bins milo-full-4m-lineage-matched-age-bins milo-full-4m-lineage-edger-parity milo-full-4m-lineage-matched-edger-parity milor-lineage-subset-parity milor-lineage-matched-subset-parity clean

setup:
	$(PYTHON) -m pip install -e ".[dev]"

test:
	$(PYTHON) -m unittest discover -s tests -p "test_*.py"

download-gateway:
	$(PYTHON) scripts/data.py download --config configs/gateway.yaml

download-info:
	$(PYTHON) scripts/data.py download --config configs/gateway.yaml --dry-run

download-gse184117:
	mkdir -p data/external
	curl -L -C - -o data/external/GSE184117_RAW.tar https://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184117/suppl/GSE184117_RAW.tar
	curl -L -C - -o data/external/GSE184117_series_matrix.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184117/matrix/GSE184117_series_matrix.txt.gz

inspect-gse184117:
	$(PYTHON) scripts/external.py inspect-archive --archive data/external/GSE184117_RAW.tar --dataset-id oliva_2022

external-gse184117-modules:
	$(PYTHON) scripts/external.py score-modules --archive data/external/GSE184117_RAW.tar --metadata data/external/GSE184117_series_matrix.txt.gz --dataset-id oliva_2022

external-gse184117-markers:
	$(PYTHON) scripts/external.py score-markers --archive data/external/GSE184117_RAW.tar --metadata data/external/GSE184117_series_matrix.txt.gz --dataset-id oliva_2022

external-gse184117-mapped:
	$(PYTHON) scripts/external.py build-mapped --archive data/external/GSE184117_RAW.tar --metadata data/external/GSE184117_series_matrix.txt.gz --dataset-id oliva_2022

external-gse184117-status:
	$(PYTHON) scripts/external.py gse184117-status --external-config configs/external_datasets.yaml

external-scanvi-reference-map:
	$(PYTHON) scripts/run_scanvi_reference_mapping.py --reference-h5ad data/processed/gateway_scvi_stratified_250k.h5ad --query-h5ad data/processed/gse184117_marker_mapped.h5ad --model-dir results/models/gateway_scanvi_reference --query-out data/processed/gse184117_scanvi_mapped.h5ad --qc-out results/tables/external_scanvi_mapping_qc.tsv --donor-features-out data/processed/gse184117_scanvi_donor_features.tsv --metadata-out results/tables/gateway_scanvi_reference_metadata.tsv --overwrite-model

external-scanvi-feature-concordance:
	$(PYTHON) scripts/external.py compare-mapped --external-config configs/external_datasets.yaml --gateway-config configs/gateway.yaml --mapped-features data/processed/gse184117_scanvi_donor_features.tsv --out results/tables/external_scanvi_feature_concordance.tsv --direct-feature-map

external-mapped-feature-concordance:
	$(PYTHON) scripts/external.py compare-mapped --external-config configs/external_datasets.yaml --gateway-config configs/gateway.yaml

external-marker-age-concordance:
	$(PYTHON) scripts/external.py compare-marker-age --gateway-config configs/gateway.yaml

inspect:
	$(PYTHON) scripts/data.py inspect --config configs/gateway.yaml

cohort:
	$(PYTHON) scripts/data.py manifest --config configs/gateway.yaml

aggregate:
	$(PYTHON) scripts/data.py aggregate --config configs/gateway.yaml

features:
	$(PYTHON) scripts/data.py features --config configs/gateway.yaml

features-augmented:
	$(PYTHON) scripts/data.py features --config configs/gateway.yaml --include-modules

age-associations:
	$(PYTHON) scripts/modeling.py age-associations --features data/processed/donor_cell_state_features.tsv --manifest data/processed/cohort_manifest.tsv --config configs/gateway.yaml

compositional-age-model:
	$(PYTHON) scripts/run_compositional_age_model.py --config configs/gateway.yaml

negative-controls:
	$(PYTHON) scripts/run_negative_controls.py --gateway-config configs/gateway.yaml --model-config configs/models.yaml

feature-family-ablation: scvi-hybrid-features
	$(PYTHON) scripts/run_feature_family_ablation.py --features data/processed/ora_scvi_hybrid_feature_matrix.tsv --manifest data/processed/cohort_manifest.tsv --model-config configs/models.yaml --gateway-config configs/gateway.yaml --models hist_gradient_boosting xgboost catboost boosted_ensemble --repeats 10 --n-permutations 20 --permutation-repeats 1

leave-context-out: scvi-hybrid-features
	$(PYTHON) scripts/run_leave_context_out.py --features data/processed/ora_scvi_hybrid_feature_matrix.tsv --manifest data/processed/cohort_manifest.tsv --model-config configs/models.yaml --gateway-config configs/gateway.yaml --models hist_gradient_boosting xgboost catboost boosted_ensemble --repeats 10 --min-train-donors 40 --min-test-donors 5

foundation-benchmark-subsets:
	$(PYTHON) scripts/foundation_benchmark.py subsets --h5ad $(FOUNDATION_BENCHMARK_H5AD) --lineage-cells $(FOUNDATION_LINEAGE_CELLS) --epithelial-cells $(FOUNDATION_EPITHELIAL_CELLS) --allcell-cells $(FOUNDATION_ALLCELL_CELLS) --overwrite

geneformer-donor-features:
	$(FOUNDATION_PYTHON) scripts/run_geneformer_benchmark.py --h5ad $(GENEFORMER_H5AD) --max-cells $(GENEFORMER_MAX_CELLS) --batch-size $(GENEFORMER_BATCH_SIZE)

geneformer-age-model: geneformer-donor-features
	$(PYTHON) scripts/modeling.py repeated --features data/processed/geneformer_donor_features.tsv --manifest data/processed/cohort_manifest.tsv --model-config configs/models.yaml --gateway-config configs/gateway.yaml --models hist_gradient_boosting xgboost catboost boosted_ensemble --repeat-performance-out results/tables/geneformer_repeated_cv_performance.tsv --summary-out results/tables/geneformer_age_model_summary.tsv --scores-out results/tables/geneformer_scores.tsv --feature-stability-out results/tables/geneformer_feature_stability.tsv --repeats 10

foundation-model-benchmark: geneformer-age-model
	$(PYTHON) scripts/foundation_benchmark.py compare --geneformer-summary results/tables/geneformer_age_model_summary.tsv --runtime results/tables/foundation_model_runtime.tsv --out results/tables/foundation_model_benchmark.tsv

aging-clock-baseline:
	$(PYTHON) scripts/run_aging_clock_baseline.py

cross-tissue-specificity:
	$(PYTHON) scripts/run_cross_tissue_specificity.py

cross-tissue-age-effects: cross-tissue-specificity
	$(PYTHON) scripts/run_cross_tissue_age_effects.py

spatial-validation-plan:
	$(PYTHON) scripts/run_spatial_validation_plan.py

perturbation-validation-plan:
	$(PYTHON) scripts/run_perturbation_validation_plan.py

perturbation-gse309325-organoid:
	$(PYTHON) scripts/run_gse309325_organoid_adapter.py

regeneration-axis-feature-map:
	$(PYTHON) scripts/run_regeneration_axis.py

regeneration-modules:
	$(PYTHON) scripts/data.py modules --config configs/gateway.yaml --gene-sets configs/regeneration_gene_sets.yaml --summary-out results/tables/regeneration_module_score_summary.tsv --coverage-out results/tables/regeneration_module_gene_coverage.tsv --donor-features-out data/processed/regeneration_donor_module_features.tsv

regeneration-module-analysis: regeneration-modules
	$(PYTHON) scripts/run_regeneration_module_analysis.py

regulatory-driver-analysis:
	$(PYTHON) scripts/run_regulatory_driver_analysis.py

niche-signaling-analysis:
	$(PYTHON) scripts/run_niche_signaling_analysis.py

regeneration-dynamics:
	$(PYTHON) scripts/run_regeneration_dynamics.py

model-ora:
	$(PYTHON) scripts/modeling.py train --features data/processed/ora_feature_matrix.tsv --manifest data/processed/cohort_manifest.tsv --config configs/models.yaml

model-ora-diagnostics:
	$(PYTHON) scripts/modeling.py diagnostics --model-config configs/models.yaml --gateway-config configs/gateway.yaml

model-ora-repeated:
	$(PYTHON) scripts/modeling.py repeated --features data/processed/ora_feature_matrix.tsv --manifest data/processed/cohort_manifest.tsv --model-config configs/models.yaml --gateway-config configs/gateway.yaml

model-ora-augmented:
	$(PYTHON) scripts/modeling.py train --features data/processed/ora_augmented_feature_matrix.tsv --manifest data/processed/cohort_manifest.tsv --config configs/models.yaml --performance-out results/tables/ora_augmented_model_performance.tsv --scores-out results/tables/augmented_donor_ora_scores.tsv --importance-out results/tables/ora_augmented_feature_importance.tsv

model-ora-candidate-repeated:
	$(PYTHON) scripts/modeling.py repeated --features data/processed/ora_augmented_feature_matrix.tsv --manifest data/processed/cohort_manifest.tsv --model-config configs/models.yaml --gateway-config configs/gateway.yaml --models hist_gradient_boosting xgboost catboost boosted_ensemble --repeat-performance-out results/tables/ora_augmented_candidate_repeated_cv_performance.tsv --summary-out results/tables/ora_augmented_candidate_repeated_cv_summary.tsv --scores-out results/tables/ora_augmented_candidate_repeated_cv_scores.tsv --feature-stability-out results/tables/ora_augmented_candidate_repeated_cv_feature_stability.tsv --repeats 10

ora-feature-set-comparison:
	$(PYTHON) scripts/modeling.py feature-set-comparison

ora-permutation-null:
	$(PYTHON) scripts/modeling.py permutation-null --observed-summary results/tables/ora_augmented_candidate_repeated_cv_summary.tsv --models hist_gradient_boosting xgboost catboost boosted_ensemble --n-permutations 50 --repeats 2

ora-nested-tuning:
	$(PYTHON) scripts/modeling.py nested-tuning --allow-fallback

ora-stacking:
	$(PYTHON) scripts/modeling.py stacking --allow-fallback

feature-interpretation:
	$(PYTHON) scripts/modeling.py interpretation --gateway-config configs/gateway.yaml

project-ndd:
	$(PYTHON) scripts/ndd.py project --gateway-config configs/gateway.yaml --model-config configs/models.yaml

project-ndd-feature-sensitivity:
	$(PYTHON) scripts/ndd.py feature-sensitivity --gateway-config configs/gateway.yaml --model-config configs/models.yaml

project-ndd-uncertainty:
	$(PYTHON) scripts/ndd.py uncertainty --gateway-config configs/gateway.yaml

project-ndd-diagnostics:
	$(PYTHON) scripts/ndd.py diagnostics --gateway-config configs/gateway.yaml

project-ndd-label-permutation:
	$(PYTHON) scripts/ndd.py label-permutation --gateway-config configs/gateway.yaml

manuscript-figures:
	$(PYTHON) scripts/build_manuscript_figures.py --tables-dir results/tables --figures-dir results/figures

publication-tables:
	$(PYTHON) scripts/reporting.py publication-tables --tables-dir results/tables --out-dir results/tables --index-out results/reports/publication_tables.md

manuscript:
	cd manuscript && if command -v latexmk >/dev/null 2>&1; then latexmk -pdf main.tex; elif command -v pdflatex >/dev/null 2>&1 && command -v bibtex >/dev/null 2>&1; then pdflatex main.tex && bibtex main && pdflatex main.tex && pdflatex main.tex; else echo "No LaTeX engine found. Install latexmk or pdflatex+bibtex to build manuscript/main.pdf."; exit 2; fi

manuscript-check:
	$(PYTHON) scripts/reporting.py manuscript-check

latent-space-audit:
	$(PYTHON) scripts/audit_latent_space.py --config configs/gateway.yaml

latent-space-recompute-plan:
	$(PYTHON) scripts/plan_latent_recompute.py --config configs/gateway.yaml

latent-space-plan: latent-space-audit

scvi-pilot:
	$(PYTHON) scripts/run_scvi_latent.py --h5ad data/raw/gateway.h5ad --out data/processed/gateway_scvi_pilot_25k.h5ad --max-cells 25000 --n-top-genes 2000 --batch-key sample_id --categorical-covariates flex_version,device_guided,sex --hvg-flavor cell_ranger --hvg-batch-key flex_version --embedding-key X_scvi --max-epochs 20 --seed 13

scvi-pilot-validation:
	$(PYTHON) scripts/validate_scvi_pilot.py --config configs/gateway.yaml

scvi-scaled-250k:
	$(PYTHON) scripts/run_scvi_latent.py --h5ad data/raw/gateway.h5ad --out data/processed/gateway_scvi_stratified_250k.h5ad --max-cells 250000 --sampling-strategy stratified --stratify-keys condition,fine_celltype,sex,flex_version,device_guided --n-top-genes 3000 --batch-key sample_id --categorical-covariates flex_version,device_guided,sex --hvg-flavor cell_ranger --hvg-batch-key flex_version --embedding-key X_scvi --max-epochs 50 --seed 13

scvi-scaled-250k-seed23:
	$(PYTHON) scripts/run_scvi_latent.py --h5ad data/raw/gateway.h5ad --out data/processed/gateway_scvi_stratified_250k_seed23.h5ad --model-dir results/models/gateway_scvi_stratified_250k_seed23 --max-cells 250000 --sampling-strategy stratified --stratify-keys condition,fine_celltype,sex,flex_version,device_guided --n-top-genes 3000 --batch-key sample_id --categorical-covariates flex_version,device_guided,sex --hvg-flavor cell_ranger --hvg-batch-key flex_version --embedding-key X_scvi --max-epochs 100 --accelerator auto --devices auto --seed 23

scvi-scaled-validation:
	$(PYTHON) scripts/validate_scvi_pilot.py --config configs/gateway.yaml --h5ad data/processed/gateway_scvi_stratified_250k.h5ad --out results/tables/scvi_scaled_250k_validation.tsv --max-validation-cells 50000 --seed 13

scvi-scaled-500k:
	$(PYTHON) scripts/run_scvi_latent.py --h5ad data/raw/gateway.h5ad --out data/processed/gateway_scvi_stratified_500k.h5ad --model-dir results/models/gateway_scvi_stratified_500k --max-cells 500000 --sampling-strategy stratified --stratify-keys condition,fine_celltype,sex,flex_version,device_guided --gene-list-file resources/scvi/full_4m_genes.txt --n-top-genes 3000 --batch-key sample_id --categorical-covariates flex_version,device_guided,sex --hvg-flavor cell_ranger --hvg-batch-key flex_version --embedding-key X_scvi --max-epochs 100 --accelerator auto --devices auto --seed 23

scvi-scaled-1m:
	$(PYTHON) scripts/run_scvi_latent.py --h5ad data/raw/gateway.h5ad --out data/processed/gateway_scvi_stratified_1m.h5ad --model-dir results/models/gateway_scvi_stratified_1m --max-cells 1000000 --sampling-strategy stratified --stratify-keys condition,fine_celltype,sex,flex_version,device_guided --gene-list-file resources/scvi/full_4m_genes.txt --n-top-genes 3000 --batch-key sample_id --categorical-covariates flex_version,device_guided,sex --hvg-flavor cell_ranger --hvg-batch-key flex_version --embedding-key X_scvi --max-epochs 100 --accelerator auto --devices auto --seed 23

scvi-reduced-4m:
	$(PYTHON) scripts/build_reduced_h5ad.py --h5ad data/raw/gateway.h5ad --gene-list-file resources/scvi/full_4m_genes.txt --out data/processed/gateway_hvg3003_4m.h5ad --chunk-dir data/processed/gateway_hvg3003_4m_chunks --chunk-size 25000 --overwrite

scvi-full-4m:
	$(PYTHON) scripts/run_scvi_latent.py --h5ad data/raw/gateway.h5ad --out data/processed/gateway_scvi_full_4m.h5ad --model-dir results/models/gateway_scvi_full_4m --gene-list-file resources/scvi/full_4m_genes.txt --n-top-genes 3000 --batch-key sample_id --categorical-covariates flex_version,device_guided,sex --hvg-flavor cell_ranger --hvg-batch-key flex_version --embedding-key X_scvi --max-epochs 100 --accelerator auto --devices auto --seed 23

scvi-full-4m-safe:
	$(PYTHON) scripts/run_scvi_latent.py --h5ad data/raw/gateway.h5ad --out data/processed/gateway_scvi_full_4m_safe.h5ad --model-dir results/models/gateway_scvi_full_4m_safe --gene-list-file resources/scvi/full_4m_genes_1500.txt --n-top-genes 1500 --batch-key sample_id --categorical-covariates flex_version,device_guided,sex --hvg-flavor cell_ranger --hvg-batch-key flex_version --embedding-key X_scvi --max-epochs 100 --accelerator auto --devices auto --seed 23

scvi-full-4m-reduced:
	$(PYTHON) scripts/run_scvi_latent.py --h5ad data/processed/gateway_hvg3003_4m.h5ad --out data/processed/gateway_scvi_full_4m_reduced.h5ad --model-dir results/models/gateway_scvi_full_4m_reduced --n-top-genes 3000 --batch-key sample_id --categorical-covariates flex_version,device_guided,sex --hvg-flavor cell_ranger --hvg-batch-key flex_version --embedding-key X_scvi --max-epochs 100 --accelerator auto --devices auto --seed 23

scvi-full-validation:
	$(PYTHON) scripts/validate_scvi_pilot.py --config configs/gateway.yaml --h5ad data/processed/gateway_scvi_full_4m_reduced.h5ad --out results/tables/scvi_full_4m_reduced_validation.tsv --max-validation-cells 100000 --seed 23

scvi-scaled-comparison:
	$(PYTHON) scripts/summarize_scvi_validations.py --validation full_4m_reduced results/tables/scvi_full_4m_reduced_validation.tsv --validation stratified_250k_seed13 results/tables/scvi_scaled_250k_validation.tsv --validation stratified_250k_seed23 results/tables/scvi_scaled_250k_seed23_validation.tsv --validation lineage_basal_neural_100k results/tables/scvi_lineage_basal_neural_validation.tsv --out results/tables/scvi_latent_validation_comparison.tsv

scvi-embedding-claim-gates:
	$(PYTHON) scripts/compare_scvi_embeddings.py --validation full_4m_reduced results/tables/scvi_full_4m_reduced_validation.tsv --validation stratified_250k_seed13 results/tables/scvi_scaled_250k_validation.tsv --validation stratified_250k_seed23 results/tables/scvi_scaled_250k_seed23_validation.tsv --validation lineage_basal_neural_100k results/tables/scvi_lineage_basal_neural_validation.tsv --summary-out results/tables/scvi_embedding_claim_gates.tsv --markers-out results/tables/scvi_embedding_marker_concordance.tsv --note-out results/reports/scvi_embedding_comparison.md

scvi-donor-features:
	$(PYTHON) scripts/scvi_donor_embedding.py features --h5ad $(SCVI_DONOR_H5AD) --features-out data/processed/scvi_donor_embedding_features.tsv --qc-out results/tables/scvi_donor_embedding_feature_qc.tsv

scvi-donor-age-model: scvi-donor-features
	$(PYTHON) scripts/modeling.py repeated --features data/processed/scvi_donor_embedding_features.tsv --manifest data/processed/cohort_manifest.tsv --model-config configs/models.yaml --gateway-config configs/gateway.yaml --models $(SCVI_DONOR_MODELS) --repeat-performance-out results/tables/scvi_donor_embedding_repeated_cv_performance.tsv --summary-out results/tables/scvi_donor_embedding_age_model_summary.tsv --scores-out results/tables/scvi_donor_embedding_scores.tsv --feature-stability-out results/tables/scvi_donor_embedding_feature_stability.tsv --repeats 10
	$(PYTHON) scripts/scvi_donor_embedding.py state-importance --feature-stability results/tables/scvi_donor_embedding_feature_stability.tsv --out results/tables/scvi_donor_embedding_state_importance.tsv

scvi-donor-comparison: scvi-donor-age-model
	$(PYTHON) scripts/scvi_donor_embedding.py compare --summary composition results/tables/ora_repeated_cv_summary.tsv --summary composition_plus_modules results/tables/ora_augmented_candidate_repeated_cv_summary.tsv --summary scvi_donor_embedding results/tables/scvi_donor_embedding_age_model_summary.tsv --out results/tables/scvi_donor_embedding_model_comparison.tsv

scvi-donor-embedding-baseline: scvi-donor-comparison

scvi-hybrid-features: features-augmented scvi-donor-features
	$(PYTHON) scripts/scvi_donor_embedding.py hybrid-features --ora-features data/processed/ora_augmented_feature_matrix.tsv --scvi-features data/processed/scvi_donor_embedding_features.tsv --out data/processed/ora_scvi_hybrid_feature_matrix.tsv

scvi-hybrid-age-model: scvi-hybrid-features
	$(PYTHON) scripts/modeling.py repeated --features data/processed/ora_scvi_hybrid_feature_matrix.tsv --manifest data/processed/cohort_manifest.tsv --model-config configs/models.yaml --gateway-config configs/gateway.yaml --models $(SCVI_DONOR_MODELS) --repeat-performance-out results/tables/ora_scvi_hybrid_repeated_cv_performance.tsv --summary-out results/tables/ora_scvi_hybrid_age_model_summary.tsv --scores-out results/tables/ora_scvi_hybrid_scores.tsv --feature-stability-out results/tables/ora_scvi_hybrid_feature_stability.tsv --repeats 10
	$(PYTHON) scripts/scvi_donor_embedding.py state-importance --feature-stability results/tables/ora_scvi_hybrid_feature_stability.tsv --out results/tables/ora_scvi_hybrid_scvi_state_importance.tsv
	$(PYTHON) scripts/scvi_donor_embedding.py family-importance --feature-stability results/tables/ora_scvi_hybrid_feature_stability.tsv --out results/tables/ora_scvi_hybrid_feature_family_importance.tsv

scvi-hybrid-comparison: scvi-hybrid-age-model
	$(PYTHON) scripts/scvi_donor_embedding.py compare --summary composition results/tables/ora_repeated_cv_summary.tsv --summary composition_plus_modules results/tables/ora_augmented_candidate_repeated_cv_summary.tsv --summary scvi_donor_embedding results/tables/scvi_donor_embedding_age_model_summary.tsv --summary ora_scvi_hybrid results/tables/ora_scvi_hybrid_age_model_summary.tsv --out results/tables/ora_scvi_hybrid_model_comparison.tsv

scvi-hybrid-benchmark: scvi-hybrid-comparison

scvi-lineage-basal-neural:
	$(PYTHON) scripts/run_scvi_latent.py --h5ad data/raw/gateway.h5ad --out data/processed/gateway_scvi_lineage_basal_neural_100k.h5ad --max-cells 100000 --sampling-strategy stratified --stratify-keys condition,fine_celltype,sex,flex_version,device_guided --include-fine-celltype-regex "basal|HBC|GBC|INP|OSN|neuro|globose|horizontal|microvillar|secretory|sustentacular" --n-top-genes 3000 --batch-key sample_id --categorical-covariates flex_version,device_guided,sex --hvg-flavor cell_ranger --hvg-batch-key flex_version --embedding-key X_scvi --max-epochs 50 --seed 17

scvi-lineage-validation:
	$(PYTHON) scripts/validate_scvi_pilot.py --config configs/gateway.yaml --h5ad data/processed/gateway_scvi_lineage_basal_neural_100k.h5ad --out results/tables/scvi_lineage_basal_neural_validation.tsv --max-validation-cells 50000 --seed 17

modules:
	$(PYTHON) scripts/data.py modules --config configs/gateway.yaml --gene-sets configs/gene_sets.yaml

published-gene-modules:
	$(PYTHON) scripts/data.py modules --config configs/gateway.yaml --gene-sets configs/external_datasets.yaml --summary-out results/tables/published_module_score_summary.tsv --coverage-out results/tables/published_module_gene_list_coverage.tsv --donor-features-out data/processed/published_donor_module_features.tsv

external-validation:
	$(PYTHON) scripts/external.py summarize-validation --external-config configs/external_datasets.yaml --gateway-config configs/gateway.yaml

external-candidate-matrix:
	$(PYTHON) scripts/external.py candidate-matrix --external-config configs/external_datasets.yaml

external-public-data-exhaustion:
	$(PYTHON) scripts/external.py public-data-exhaustion --external-config configs/external_datasets.yaml

external-evidence:
	$(PYTHON) scripts/external.py evidence --external-config configs/external_datasets.yaml

external-feature-harmonization:
	$(PYTHON) scripts/external.py validate-feature-matrix --external-config configs/external_datasets.yaml --gateway-config configs/gateway.yaml

pseudobulk:
	$(PYTHON) scripts/pseudobulk.py targeted --config configs/gateway.yaml --gene-sets configs/gene_sets.yaml

pseudobulk-genomewide:
	$(PYTHON) scripts/pseudobulk.py export-genomewide --config configs/gateway.yaml

pseudobulk-genomewide-qc:
	$(PYTHON) scripts/pseudobulk.py qc --config configs/gateway.yaml

pseudobulk-genomewide-edger:
	$(RSCRIPT) scripts/run_genomewide_edger.R --counts data/processed/pseudobulk_genomewide_counts.tsv.gz --metadata data/processed/pseudobulk_genomewide_metadata.tsv --manifest data/processed/cohort_manifest.tsv --out results/tables/pseudobulk_genomewide_edger.tsv.gz --summary-out results/tables/pseudobulk_genomewide_edger_summary.tsv

pseudobulk-genomewide-edger-matched:
	$(RSCRIPT) scripts/run_genomewide_edger.R --counts data/processed/pseudobulk_genomewide_counts.tsv.gz --metadata data/processed/pseudobulk_genomewide_metadata.tsv --manifest data/processed/cohort_manifest.tsv --chemistry flex_v2 --collection-method device --out results/tables/pseudobulk_genomewide_edger_matched_flex_v2_device.tsv.gz --summary-out results/tables/pseudobulk_genomewide_edger_matched_flex_v2_device_summary.tsv

pseudobulk-genomewide-limma:
	$(RSCRIPT) scripts/run_genomewide_limma_voom.R --counts data/processed/pseudobulk_genomewide_counts.tsv.gz --metadata data/processed/pseudobulk_genomewide_metadata.tsv --manifest data/processed/cohort_manifest.tsv --out results/tables/pseudobulk_genomewide_limma_voom.tsv.gz --summary-out results/tables/pseudobulk_genomewide_limma_voom_summary.tsv

pseudobulk-genomewide-limma-matched:
	$(RSCRIPT) scripts/run_genomewide_limma_voom.R --counts data/processed/pseudobulk_genomewide_counts.tsv.gz --metadata data/processed/pseudobulk_genomewide_metadata.tsv --manifest data/processed/cohort_manifest.tsv --chemistry flex_v2 --collection-method device --out results/tables/pseudobulk_genomewide_limma_voom_matched_flex_v2_device.tsv.gz --summary-out results/tables/pseudobulk_genomewide_limma_voom_matched_flex_v2_device_summary.tsv

pseudobulk-genomewide-de-summary:
	$(PYTHON) scripts/pseudobulk.py de-summary --config configs/gateway.yaml

pseudobulk-genomewide-de-summary-matched:
	$(PYTHON) scripts/pseudobulk.py de-summary --config configs/gateway.yaml --de results/tables/pseudobulk_genomewide_edger_matched_flex_v2_device.tsv.gz --run-summary results/tables/pseudobulk_genomewide_edger_matched_flex_v2_device_summary.tsv --summary-out results/tables/pseudobulk_genomewide_de_summary_matched_flex_v2_device.tsv --top-hits-out results/tables/pseudobulk_genomewide_de_top_hits_matched_flex_v2_device.tsv

pseudobulk-genomewide-limma-de-summary:
	$(PYTHON) scripts/pseudobulk.py de-summary --config configs/gateway.yaml --de results/tables/pseudobulk_genomewide_limma_voom.tsv.gz --run-summary results/tables/pseudobulk_genomewide_limma_voom_summary.tsv --summary-out results/tables/pseudobulk_genomewide_limma_voom_de_summary.tsv --top-hits-out results/tables/pseudobulk_genomewide_limma_voom_de_top_hits.tsv

pseudobulk-genomewide-limma-de-summary-matched:
	$(PYTHON) scripts/pseudobulk.py de-summary --config configs/gateway.yaml --de results/tables/pseudobulk_genomewide_limma_voom_matched_flex_v2_device.tsv.gz --run-summary results/tables/pseudobulk_genomewide_limma_voom_matched_flex_v2_device_summary.tsv --summary-out results/tables/pseudobulk_genomewide_limma_voom_de_summary_matched_flex_v2_device.tsv --top-hits-out results/tables/pseudobulk_genomewide_limma_voom_de_top_hits_matched_flex_v2_device.tsv

pseudobulk-genomewide-de-audit:
	$(PYTHON) scripts/pseudobulk.py de-audit --config configs/gateway.yaml

pseudobulk-genomewide-de-audit-matched:
	$(PYTHON) scripts/pseudobulk.py de-audit --config configs/gateway.yaml --de results/tables/pseudobulk_genomewide_edger_matched_flex_v2_device.tsv.gz --run-summary results/tables/pseudobulk_genomewide_edger_matched_flex_v2_device_summary.tsv --audit-out results/tables/pseudobulk_genomewide_de_audit_matched_flex_v2_device.tsv --donor-balance-out results/tables/pseudobulk_genomewide_donor_balance_matched_flex_v2_device.tsv --matched-feasibility-out results/tables/pseudobulk_genomewide_matched_feasibility_matched_flex_v2_device.tsv

pseudobulk-genomewide-limma-de-audit:
	$(PYTHON) scripts/pseudobulk.py de-audit --config configs/gateway.yaml --de results/tables/pseudobulk_genomewide_limma_voom.tsv.gz --run-summary results/tables/pseudobulk_genomewide_limma_voom_summary.tsv --audit-out results/tables/pseudobulk_genomewide_limma_voom_de_audit.tsv --donor-balance-out results/tables/pseudobulk_genomewide_limma_voom_donor_balance.tsv --matched-feasibility-out results/tables/pseudobulk_genomewide_limma_voom_matched_feasibility.tsv

pseudobulk-genomewide-limma-de-audit-matched:
	$(PYTHON) scripts/pseudobulk.py de-audit --config configs/gateway.yaml --de results/tables/pseudobulk_genomewide_limma_voom_matched_flex_v2_device.tsv.gz --run-summary results/tables/pseudobulk_genomewide_limma_voom_matched_flex_v2_device_summary.tsv --audit-out results/tables/pseudobulk_genomewide_limma_voom_de_audit_matched_flex_v2_device.tsv --donor-balance-out results/tables/pseudobulk_genomewide_limma_voom_donor_balance_matched_flex_v2_device.tsv --matched-feasibility-out results/tables/pseudobulk_genomewide_limma_voom_matched_feasibility_matched_flex_v2_device.tsv

pseudobulk-covariate-de:
	$(PYTHON) scripts/pseudobulk.py covariate-de --config configs/gateway.yaml

ora-sensitivity:
	$(PYTHON) scripts/modeling.py sensitivity --gateway-config configs/gateway.yaml --model-config configs/models.yaml

ora-sensitivity-rf:
	$(PYTHON) scripts/modeling.py sensitivity --gateway-config configs/gateway.yaml --model-config configs/models.yaml --models random_forest

model-card:
	$(PYTHON) scripts/reporting.py model-card --gateway-config configs/gateway.yaml

output-provenance:
	$(PYTHON) scripts/reporting.py output-provenance --gateway-config configs/gateway.yaml --command-manifest configs/command_manifest.yaml

release-manifest:
	$(PYTHON) scripts/build_release_manifest.py --config configs/gateway.yaml --command-manifest configs/command_manifest.yaml

environment-report:
	$(PYTHON) scripts/build_environment_report.py --config configs/gateway.yaml

archive-review-package: release-manifest
	$(PYTHON) scripts/build_archive_review_package.py --config configs/gateway.yaml

local-light:
	$(PYTHON) -m ruff check .
	$(MAKE) test
	$(MAKE) external-validation
	$(MAKE) external-candidate-matrix
	$(MAKE) external-public-data-exhaustion
	$(MAKE) external-gse184117-status
	$(MAKE) publication-tables
	$(MAKE) manuscript-figures
	$(MAKE) manuscript-check
	$(MAKE) output-provenance

local-medium:
	$(MAKE) cohort
	$(MAKE) aggregate
	$(MAKE) features
	$(MAKE) features-augmented
	$(MAKE) age-associations
	$(MAKE) compositional-age-model
	$(MAKE) negative-controls
	$(MAKE) model-ora-diagnostics
	$(MAKE) model-ora-repeated
	$(MAKE) model-ora-candidate-repeated
	$(MAKE) ora-feature-set-comparison
	$(MAKE) feature-family-ablation
	$(MAKE) leave-context-out
	$(MAKE) ora-permutation-null
	$(MAKE) ora-nested-tuning
	$(MAKE) ora-stacking
	$(MAKE) feature-interpretation
	$(MAKE) external-evidence
	$(MAKE) project-ndd-feature-sensitivity
	$(MAKE) project-ndd-uncertainty
	$(MAKE) project-ndd-diagnostics
	$(MAKE) project-ndd-label-permutation
	$(MAKE) local-light

remote-heavy:
	$(MAKE) scvi-reduced-4m
	$(MAKE) scvi-full-4m-reduced
	$(MAKE) scvi-full-validation
	$(MAKE) scvi-scaled-250k-seed23
	$(MAKE) scvi-scaled-validation
	$(MAKE) scvi-scaled-comparison
	$(MAKE) scvi-embedding-claim-gates
	$(MAKE) scvi-donor-embedding-baseline
	$(MAKE) scvi-hybrid-benchmark
	$(MAKE) milo-full-4m-lineage
	$(MAKE) milo-full-4m-lineage-matched
	$(MAKE) milo-full-4m-lineage-programs
	$(MAKE) milo-full-4m-lineage-matched-programs
	$(MAKE) milo-full-4m-lineage-age-bins
	$(MAKE) milo-full-4m-lineage-matched-age-bins
	$(MAKE) milo-full-4m-lineage-edger-parity
	$(MAKE) milo-full-4m-lineage-matched-edger-parity
	$(MAKE) milor-lineage-subset-parity
	$(MAKE) milor-lineage-matched-subset-parity

submission-freeze:
	$(MAKE) environment-report
	$(MAKE) release-manifest
	$(MAKE) archive-review-package
	$(MAKE) external-gse184117-status
	$(MAKE) external-public-data-exhaustion
	$(MAKE) publication-tables
	$(MAKE) manuscript-figures
	$(MAKE) manuscript-check
	$(MAKE) output-provenance

all-summary: external-validation external-candidate-matrix external-gse184117-modules external-gse184117-markers external-gse184117-mapped external-marker-age-concordance external-mapped-feature-concordance external-evidence model-ora-diagnostics ora-feature-set-comparison feature-family-ablation leave-context-out ora-permutation-null ora-nested-tuning ora-stacking feature-interpretation scvi-donor-comparison scvi-hybrid-comparison pseudobulk-genomewide-de-summary pseudobulk-genomewide-de-audit pseudobulk-genomewide-de-summary-matched pseudobulk-genomewide-de-audit-matched pseudobulk-genomewide-limma-de-summary pseudobulk-genomewide-limma-de-audit pseudobulk-genomewide-limma-de-summary-matched pseudobulk-genomewide-limma-de-audit-matched project-ndd-feature-sensitivity project-ndd-uncertainty project-ndd-diagnostics project-ndd-label-permutation model-card latent-space-audit latent-space-recompute-plan output-provenance

remote-full-scvi:
	tools/remote_full_scvi.sh $(REMOTE_FULL_SCVI_ACTION)

milo:
	$(PYTHON) scripts/run_milo_pilot.py --h5ad $(MILO_H5AD) --manifest data/processed/cohort_manifest.tsv --n-neighborhoods $(MILO_NEIGHBORHOODS) --n-neighbors $(MILO_NEIGHBORS) --seed-stratify-columns coarse_celltype,fine_celltype --out results/tables/milo_pilot_neighborhood_da.tsv --summary-out results/tables/milo_pilot_summary.tsv

milo-lineage:
	$(PYTHON) scripts/run_milo_pilot.py --h5ad $(MILO_LINEAGE_H5AD) --manifest data/processed/cohort_manifest.tsv --include-coarse-regex "Resp_HBC|Olf_INPs|Olf_iOSNs|Olf_mOSNs|Olf_Sus" --n-neighborhoods $(MILO_NEIGHBORHOODS) --n-neighbors $(MILO_NEIGHBORS) --seed-stratify-columns coarse_celltype,fine_celltype --out results/tables/milo_lineage_neighborhood_da.tsv --summary-out results/tables/milo_lineage_summary.tsv

milo-secretory:
	$(PYTHON) scripts/run_milo_pilot.py --h5ad $(MILO_H5AD) --manifest data/processed/cohort_manifest.tsv --include-coarse-regex "Resp_Secretory|Olf_Sus|Bowman_Gland" --n-neighborhoods $(MILO_NEIGHBORHOODS) --n-neighbors $(MILO_NEIGHBORS) --seed-stratify-columns coarse_celltype,fine_celltype --out results/tables/milo_secretory_neighborhood_da.tsv --summary-out results/tables/milo_secretory_summary.tsv

milo-full-4m:
	$(PYTHON) scripts/run_milo_pilot.py --h5ad $(MILO_FULL_H5AD) --manifest data/processed/cohort_manifest.tsv --n-neighborhoods $(MILO_FULL_NEIGHBORHOODS) --n-neighbors $(MILO_FULL_NEIGHBORS) --min-donors 30 --seed-stratify-columns coarse_celltype,fine_celltype --out results/tables/milo_full_4m_neighborhood_da.tsv --summary-out results/tables/milo_full_4m_summary.tsv

milo-full-4m-lineage:
	$(PYTHON) scripts/run_milo_pilot.py --h5ad $(MILO_FULL_H5AD) --manifest data/processed/cohort_manifest.tsv --include-coarse-regex "Resp_HBC|Olf_INPs|Olf_iOSNs|Olf_mOSNs|Olf_Sus" --n-neighborhoods $(MILO_FULL_NEIGHBORHOODS) --n-neighbors $(MILO_FULL_NEIGHBORS) --min-donors 30 --seed-stratify-columns coarse_celltype,fine_celltype --out results/tables/milo_full_4m_lineage_neighborhood_da.tsv --summary-out results/tables/milo_full_4m_lineage_summary.tsv

milo-full-4m-secretory:
	$(PYTHON) scripts/run_milo_pilot.py --h5ad $(MILO_FULL_H5AD) --manifest data/processed/cohort_manifest.tsv --include-coarse-regex "Resp_Secretory|Olf_Sus|Bowman_Gland" --n-neighborhoods $(MILO_FULL_NEIGHBORHOODS) --n-neighbors $(MILO_FULL_NEIGHBORS) --min-donors 30 --seed-stratify-columns coarse_celltype,fine_celltype --out results/tables/milo_full_4m_secretory_neighborhood_da.tsv --summary-out results/tables/milo_full_4m_secretory_summary.tsv

milo-full-4m-matched:
	$(PYTHON) scripts/run_milo_pilot.py --h5ad $(MILO_FULL_H5AD) --manifest data/processed/cohort_manifest.tsv --donor-query "$(MILO_MATCHED_DONOR_QUERY)" --n-neighborhoods $(MILO_FULL_NEIGHBORHOODS) --n-neighbors $(MILO_FULL_NEIGHBORS) --min-donors $(MILO_MATCHED_MIN_DONORS) --seed-stratify-columns coarse_celltype,fine_celltype --out results/tables/milo_full_4m_matched_neighborhood_da.tsv --summary-out results/tables/milo_full_4m_matched_summary.tsv

milo-full-4m-lineage-matched:
	$(PYTHON) scripts/run_milo_pilot.py --h5ad $(MILO_FULL_H5AD) --manifest data/processed/cohort_manifest.tsv --donor-query "$(MILO_MATCHED_DONOR_QUERY)" --include-coarse-regex "Resp_HBC|Olf_INPs|Olf_iOSNs|Olf_mOSNs|Olf_Sus" --n-neighborhoods $(MILO_FULL_NEIGHBORHOODS) --n-neighbors $(MILO_FULL_NEIGHBORS) --min-donors $(MILO_MATCHED_MIN_DONORS) --seed-stratify-columns coarse_celltype,fine_celltype --out results/tables/milo_full_4m_lineage_matched_neighborhood_da.tsv --summary-out results/tables/milo_full_4m_lineage_matched_summary.tsv

milo-full-4m-secretory-matched:
	$(PYTHON) scripts/run_milo_pilot.py --h5ad $(MILO_FULL_H5AD) --manifest data/processed/cohort_manifest.tsv --donor-query "$(MILO_MATCHED_DONOR_QUERY)" --include-coarse-regex "Resp_Secretory|Olf_Sus|Bowman_Gland" --n-neighborhoods $(MILO_FULL_NEIGHBORHOODS) --n-neighbors $(MILO_FULL_NEIGHBORS) --min-donors $(MILO_MATCHED_MIN_DONORS) --seed-stratify-columns coarse_celltype,fine_celltype --out results/tables/milo_full_4m_secretory_matched_neighborhood_da.tsv --summary-out results/tables/milo_full_4m_secretory_matched_summary.tsv

milo-full-4m-lineage-programs:
	$(PYTHON) scripts/score_milo_neighborhood_programs.py --h5ad $(MILO_FULL_H5AD) --memberships results/tables/milo_full_4m_lineage_memberships.tsv --da-table results/tables/milo_full_4m_lineage_neighborhood_da.tsv --run-name lineage_full --scores-out results/tables/milo_full_4m_lineage_program_scores.tsv --summary-out results/tables/milo_full_4m_lineage_program_summary.tsv --coverage-out results/tables/milo_full_4m_lineage_program_coverage.tsv

milo-full-4m-lineage-matched-programs:
	$(PYTHON) scripts/score_milo_neighborhood_programs.py --h5ad $(MILO_FULL_H5AD) --memberships results/tables/milo_full_4m_lineage_matched_memberships.tsv --da-table results/tables/milo_full_4m_lineage_matched_neighborhood_da.tsv --run-name lineage_matched --scores-out results/tables/milo_full_4m_lineage_matched_program_scores.tsv --summary-out results/tables/milo_full_4m_lineage_matched_program_summary.tsv --coverage-out results/tables/milo_full_4m_lineage_matched_program_coverage.tsv

milo-full-4m-lineage-age-bins:
	$(PYTHON) scripts/summarize_milo_age_bins.py --memberships results/tables/milo_full_4m_lineage_memberships.tsv --manifest data/processed/cohort_manifest.tsv --da-table results/tables/milo_full_4m_lineage_neighborhood_da.tsv --run-name lineage_full --neighborhoods-out results/tables/milo_full_4m_lineage_age_bin_neighborhoods.tsv --summary-out results/tables/milo_full_4m_lineage_age_bin_summary.tsv

milo-full-4m-lineage-matched-age-bins:
	$(PYTHON) scripts/summarize_milo_age_bins.py --memberships results/tables/milo_full_4m_lineage_matched_memberships.tsv --manifest data/processed/cohort_manifest.tsv --donor-query "$(MILO_MATCHED_DONOR_QUERY)" --da-table results/tables/milo_full_4m_lineage_matched_neighborhood_da.tsv --run-name lineage_matched --neighborhoods-out results/tables/milo_full_4m_lineage_matched_age_bin_neighborhoods.tsv --summary-out results/tables/milo_full_4m_lineage_matched_age_bin_summary.tsv

milo-full-4m-lineage-edger-parity:
	$(PYTHON) scripts/export_milo_neighborhood_counts.py --memberships results/tables/milo_full_4m_lineage_memberships.tsv --manifest data/processed/cohort_manifest.tsv --counts-out data/processed/milo_full_4m_lineage_neighborhood_counts.tsv --design-out data/processed/milo_full_4m_lineage_neighborhood_design.tsv --summary-out results/tables/milo_full_4m_lineage_count_export_summary.tsv
	$(RSCRIPT) scripts/run_milo_edger_parity.R --counts data/processed/milo_full_4m_lineage_neighborhood_counts.tsv --design data/processed/milo_full_4m_lineage_neighborhood_design.tsv --out results/tables/milo_full_4m_lineage_edger_parity.tsv --summary-out results/tables/milo_full_4m_lineage_edger_parity_run_summary.tsv
	$(PYTHON) scripts/summarize_milo_edger_parity.py --python-da results/tables/milo_full_4m_lineage_neighborhood_da.tsv --edger-da results/tables/milo_full_4m_lineage_edger_parity.tsv --run-name lineage_full --comparison-out results/tables/milo_full_4m_lineage_edger_parity_comparison.tsv --summary-out results/tables/milo_full_4m_lineage_edger_parity_summary.tsv

milo-full-4m-lineage-matched-edger-parity:
	$(PYTHON) scripts/export_milo_neighborhood_counts.py --memberships results/tables/milo_full_4m_lineage_matched_memberships.tsv --manifest data/processed/cohort_manifest.tsv --donor-query "$(MILO_MATCHED_DONOR_QUERY)" --counts-out data/processed/milo_full_4m_lineage_matched_neighborhood_counts.tsv --design-out data/processed/milo_full_4m_lineage_matched_neighborhood_design.tsv --summary-out results/tables/milo_full_4m_lineage_matched_count_export_summary.tsv
	$(RSCRIPT) scripts/run_milo_edger_parity.R --counts data/processed/milo_full_4m_lineage_matched_neighborhood_counts.tsv --design data/processed/milo_full_4m_lineage_matched_neighborhood_design.tsv --out results/tables/milo_full_4m_lineage_matched_edger_parity.tsv --summary-out results/tables/milo_full_4m_lineage_matched_edger_parity_run_summary.tsv
	$(PYTHON) scripts/summarize_milo_edger_parity.py --python-da results/tables/milo_full_4m_lineage_matched_neighborhood_da.tsv --edger-da results/tables/milo_full_4m_lineage_matched_edger_parity.tsv --run-name lineage_matched --comparison-out results/tables/milo_full_4m_lineage_matched_edger_parity_comparison.tsv --summary-out results/tables/milo_full_4m_lineage_matched_edger_parity_summary.tsv

milor-lineage-subset-parity:
	$(PYTHON) scripts/export_milor_subset.py --h5ad $(MILO_FULL_H5AD) --manifest data/processed/cohort_manifest.tsv --max-cells $(MILOR_SUBSET_CELLS) --metadata-out data/processed/milor_lineage_subset_metadata.tsv --embedding-out data/processed/milor_lineage_subset_embedding.tsv --summary-out results/tables/milor_lineage_subset_export_summary.tsv
	$(MILOR_RSCRIPT) scripts/run_milor_parity.R --metadata data/processed/milor_lineage_subset_metadata.tsv --embedding data/processed/milor_lineage_subset_embedding.tsv --run-name lineage_full_subset --out results/tables/milor_lineage_subset_da.tsv --summary-out results/tables/milor_lineage_subset_summary.tsv

milor-lineage-matched-subset-parity:
	$(PYTHON) scripts/export_milor_subset.py --h5ad $(MILO_FULL_H5AD) --manifest data/processed/cohort_manifest.tsv --donor-query "$(MILO_MATCHED_DONOR_QUERY)" --max-cells $(MILOR_MATCHED_SUBSET_CELLS) --metadata-out data/processed/milor_lineage_matched_subset_metadata.tsv --embedding-out data/processed/milor_lineage_matched_subset_embedding.tsv --summary-out results/tables/milor_lineage_matched_subset_export_summary.tsv
	$(MILOR_RSCRIPT) scripts/run_milor_parity.R --metadata data/processed/milor_lineage_matched_subset_metadata.tsv --embedding data/processed/milor_lineage_matched_subset_embedding.tsv --run-name lineage_matched_subset --out results/tables/milor_lineage_matched_subset_da.tsv --summary-out results/tables/milor_lineage_matched_subset_summary.tsv

clean:
	rm -rf .pytest_cache .ruff_cache build dist *.egg-info
