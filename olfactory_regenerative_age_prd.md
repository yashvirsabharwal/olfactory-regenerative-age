# PRD and Methodology: Olfactory Regenerative Age (ORA)

## Working paper title

**Olfactory Regenerative Age: a population-scale single-cell model of adult human olfactory neurogenesis and aging**

Alternative titles:

1. **Population-scale mapping of adult human olfactory neurogenesis reveals age-associated regenerative bottlenecks**
2. **A single-cell regenerative aging clock for the human olfactory neuroepithelium**
3. **Living human olfactory neurons reveal a regenerative aging axis linked to impaired neurogenesis and neurodegenerative disease programs**

## One-sentence paper claim

A large, newly released human olfactory epithelium single-cell atlas enables a donor-level, interpretable machine-learning model of olfactory regenerative aging that identifies where the adult human olfactory neurogenic lineage becomes bottlenecked with age and links those shifts to inflammatory, stem-cell dysfunction, neuronal maturation, and neurodegenerative-disease-relevant programs.

## Core concept

The Gateway atlas contains millions of cells from hundreds of donors, with healthy controls spanning a broad adult age range and fine-resolution labels across basal cells, neuronal progenitors, immature neurons, mature neurons, stressed mature neurons, sustentacular cells, respiratory epithelium, glandular cells, and immune populations. This creates a rare opportunity to move beyond descriptive single-cell aging analysis and build a **quantitative model of regenerative capacity** in a living human neuroepithelium.

The central output is a donor-level score:

**ORA = Olfactory Regenerative Age**

ORA is an interpretable predicted-age score learned from cell-state composition, gene-program activity, and lineage-position features in healthy donors. The residual after adjusting for chronological age and technical covariates is:

**ORAA = Olfactory Regenerative Age Acceleration**

A positive ORAA indicates an olfactory neuroepithelium that looks older than expected; a negative ORAA indicates a more youthful regenerative profile.

A second biological score is derived from the model:

**NCI = Neurogenic Capacity Index**

NCI summarizes the abundance and transcriptional state of the regenerative lineage from horizontal basal cells to intermediate neuronal progenitors, immature olfactory sensory neurons, mature olfactory sensory neurons, and stressed mature neurons.

## Why this is a fast and publishable direction

This project avoids the main weakness of the Gateway preprint's disease analysis: the disease cohort is small. Instead, the primary discovery cohort is the large healthy-control subset. The neurodegenerative disease donors are used only as exploratory disease anchoring, not as the basis of a diagnostic classifier.

The paper becomes novel because it is not simply re-plotting the atlas by age. It creates a reproducible, donor-level model of regenerative aging, identifies bottlenecks along the olfactory neurogenic trajectory, validates against prior olfactory aging and Alzheimer's olfactory biopsy studies, and provides an analysis framework others can reuse.

## Main hypotheses

### H1: Aging changes the composition of the olfactory neurogenic lineage

Older healthy donors will show altered proportions of basal cells, activated basal cells, neuronal progenitors, immature olfactory sensory neurons, mature olfactory sensory neurons, and stressed mature olfactory sensory neurons.

Expected directions to test, not assume:

- Lower progenitor-to-neuron output with age.
- Lower immature-to-mature neuronal maturation efficiency.
- Increased stressed mature neuronal state.
- Increased respiratory epithelial or metaplastic-like contribution.
- Increased inflammatory immune or epithelial stress signatures.

### H2: Aging creates a regenerative bottleneck, not a uniform loss across all stages

The largest age-associated shift should occur at a specific transition, for example:

- Quiescent HBC -> activated HBC
- Activated HBC -> INP
- INP -> iOSN
- iOSN -> fully mature mOSN
- Fully mature mOSN -> stressed mOSN

The paper should identify which transition is most affected and support it with abundance, pseudotime-density, and gene-program evidence.

### H3: Donor-level ORA can be predicted from interpretable single-cell features

A regularized age model trained on healthy donors should predict chronological age better than null models and should reveal a sparse set of biological features explaining olfactory aging.

### H4: ORA overlaps with presbyosmia and neurodegeneration-related olfactory signatures

Even though Gateway disease sample size is small, ORA-related gene programs should overlap with prior presbyosmia and Alzheimer's olfactory biopsy datasets, especially stem-cell inflammation, cytokine response, TP63-associated differentiation blockade, immune activation, and neuronal injury signatures.

## Expected manuscript structure

### Abstract

- Living human olfactory epithelium supports adult neurogenesis, but population-scale maps of regenerative aging are lacking.
- Use Gateway 4M atlas to model age-associated changes in healthy human olfactory neuroepithelium.
- Define ORA and NCI from cell-state composition, pathway activity, and lineage-position features.
- Identify a regenerative bottleneck and gene programs linked to inflammation, stem-cell dysfunction, neuronal maturation, cilia/synapse biology, proteostasis, and mitochondrial/lysosomal stress.
- Validate signatures against external olfactory aging and Alzheimer's olfactory biopsy datasets.
- Provide open pipeline and reusable donor-level feature table.

### Introduction

1. Adult human olfactory epithelium is a living neurogenic tissue.
2. Olfactory function declines with age and is relevant to neurodegeneration.
3. Prior human olfactory single-cell studies were small and could not quantify population-scale aging.
4. Gateway enables a donor-level computational model of regenerative aging.
5. The contribution: ORA model, bottleneck analysis, gene programs, validation.

### Results

1. **Gateway healthy donors enable population-scale analysis of human olfactory aging**
2. **Age-associated differential abundance reveals shifts across basal, neuronal, epithelial, and immune states**
3. **Lineage-density modeling identifies a regenerative bottleneck in adult olfactory neurogenesis**
4. **ORA predicts donor age from interpretable single-cell features**
5. **Gene programs underlying ORA implicate inflammation, stem-cell dysfunction, neuronal maturation, and proteostasis**
6. **ORA-associated programs validate in external presbyosmia and Alzheimer's olfactory biopsy datasets**
7. **Exploratory disease anchoring: neurodegenerative disease donors shift toward older ORA profiles**

### Discussion

- Aging of the human olfactory neuroepithelium is a systems-level regenerative phenotype, not just loss of mature neurons.
- ORA provides a quantitative framework for comparing donor-level regenerative capacity.
- Findings support olfactory epithelium as a living tissue model for neurogenesis, olfactory loss, and possibly early neurodegenerative biology.
- Limitations: cross-sectional, transcriptomic only, possible residual chemistry/device/site confounding, no direct smell-test phenotype in Gateway unless metadata include it, disease cohort small.

---

# Product Requirements Document

## Product name

`olfactory-regenerative-age`

## Product type

Reproducible computational biology and machine-learning research pipeline.

## Primary user

Researcher using Codex to build a reproducible analysis package from public single-cell data.

## Secondary users

- Computational biologists reproducing results.
- Reviewers evaluating statistical validity.
- Neurodegeneration or olfaction researchers reusing ORA features.
- Future users applying ORA to new olfactory biopsy datasets.

## Goal

Build a complete, reproducible pipeline that ingests the Gateway CELLxGENE atlas and external olfactory datasets, computes donor-level regenerative aging features, fits and validates ORA models, identifies cell-state and pathway drivers of aging, and outputs manuscript-ready figures and tables.

## Non-goals

- Do not build a diagnostic AD or PD classifier from Gateway alone.
- Do not claim causal aging mechanisms from cross-sectional data.
- Do not claim real lineage flux from abundance alone; use terms like "cross-sectional lineage-density proxy" unless longitudinal data are available.
- Do not require full reprocessing from FASTQ.
- Do not require full 4M-cell deep learning on a laptop.

## Success criteria

### Scientific success

The project is successful if it produces:

1. A reproducible healthy-donor feature table.
2. A cross-validated ORA model with performance above null permutations.
3. Statistically robust age-associated cell-state abundance changes after covariate adjustment.
4. A lineage-density or trajectory analysis identifying candidate bottleneck regions.
5. Interpretable gene programs associated with ORA.
6. At least one external validation using published olfactory aging, presbyosmia, or Alzheimer's olfactory biopsy data.
7. A coherent set of figures sufficient for a preprint.

### Engineering success

The project is successful if:

1. A new machine can reproduce all results from a small config file.
2. Full data outputs are not committed to Git.
3. Each pipeline step writes versioned intermediate artifacts.
4. Laptop mode works with 32 GB RAM using backed AnnData, chunked aggregation, and subsampling.
5. External-compute mode supports full-neighborhood or full-subset analyses.
6. All model outputs are donor-split, never cell-split, to avoid leakage.

---

# Data sources

## Primary dataset

### Gateway 4M human olfactory epithelium atlas

Source: CELLxGENE collection provided by the user.

Primary use:

- Discovery cohort: healthy controls only.
- Exploratory anchoring: AD and PD donors only after the healthy ORA model is frozen.

Required metadata fields to inspect:

- donor_id
- sample_id
- age
- sex
- race / ethnicity
- disease condition
- FLEX chemistry version
- collection method / device usage
- coarse cell type
- fine cell type
- QC metrics such as n_counts, n_genes, percent_mito
- site if available
- UMAP/scANVI embedding if included

Expected cell states of interest:

- Quiescent HBC
- Activated HBC
- Cycling HBC
- Early INP
- Late INP
- Early iOSN
- Late iOSN
- Early mature mOSN
- Fully mature mOSN
- Stressed mOSN
- Sustentacular cells
- Respiratory ciliated cells
- Secretory epithelial cells
- Bowman gland / duct cells
- Immune populations

## External validation datasets

### Durante et al. 2020, Nature Neuroscience

Use: lineage and marker validation for adult human olfactory neurogenesis.

Expected use cases:

- Confirm conserved ordering of HBC -> progenitor -> iOSN -> mOSN trajectory.
- Compare marker modules across datasets.
- Possibly transfer Gateway ORA lineage modules onto Durante cells if metadata permit.

### Oliva et al. 2022, Journal of Clinical Investigation

Use: presbyosmia and aging validation.

Expected use cases:

- Test whether ORA-associated HBC inflammatory/stress modules are elevated in presbyosmic basal cells.
- Compare TP63 / cytokine / differentiation-blockade signatures.
- Validate that ORA captures biologically meaningful olfactory aging rather than only chronological age.

### D'Anniballe et al. 2026, Nature Communications

Use: Alzheimer's olfactory biopsy validation.

Expected use cases:

- Test whether preclinical and clinical AD olfactory biopsy samples show elevated ORA-associated inflammatory and neuronal-injury modules.
- Assess whether AD disease-stage signatures overlap with ORA gene programs.
- Keep this as validation, not training.

### Optional external datasets

- Human olfactory datasets from Goldstein lab not already above.
- SEA-AD postmortem brain atlas for exploratory overlap of ORA genes with brain aging / AD programs.
- Human Protein Atlas / GTEx tissue specificity for contextualizing neuronal genes.

---

# Technical constraints

## Local machine

Laptop: Apple M-series, 32 GB RAM.

Laptop mode should support:

- Downloading and inspecting metadata.
- Backed AnnData reads.
- Chunked donor-level cell counts.
- Donor-level pseudobulk aggregation.
- Donor-level feature engineering.
- Age model fitting.
- Most plotting.
- Downsampled trajectory analysis.

Laptop mode should avoid:

- Loading the full 4M x 18k expression matrix into dense memory.
- Training scVI from scratch on all cells.
- Full 4M-cell kNN graph construction unless precomputed embeddings and memory-efficient tools are available.
- Full unchunked UCell scoring across all cells and all pathways.

## External compute

Use external compute for:

- Full-neighborhood Milo analysis on all relevant lineage cells.
- Full scVI/scANVI retraining if needed.
- cNMF on large cell subsets.
- Large external dataset harmonization.
- Multiple sensitivity reruns.

Recommended external compute:

- 32 to 64 CPU cores.
- 128 to 512 GB RAM.
- 1 to 2 TB local SSD scratch.
- GPU optional, useful only if retraining deep models.

---

# Recommended software stack

## Python

- Python >= 3.11
- scanpy
- anndata
- scipy
- numpy
- pandas or polars
- scikit-learn
- statsmodels
- xgboost or lightgbm
- shap
- matplotlib
- seaborn optional for internal exploration, but final figures can use matplotlib
- cellxgene-census
- decoupler or gseapy
- pyarrow
- joblib
- typer or click
- pydantic
- tqdm

## R

- R >= 4.3
- SingleCellExperiment
- miloR
- edgeR
- limma
- dreamlet / variancePartition optional
- DESeq2 optional
- scCODA can be run from Python if preferred

## Workflow

- Snakemake recommended.
- Makefile wrapper for convenience.
- YAML configs for dataset paths and model settings.

---

# Repository structure

```text
olfactory-regenerative-age/
  README.md
  LICENSE
  pyproject.toml
  environment.yml
  renv.lock                      # optional for R reproducibility
  Makefile
  .gitignore

  configs/
    gateway.yaml
    external_datasets.yaml
    features.yaml
    gene_sets.yaml
    models.yaml
    plotting.yaml

  data/
    raw/                         # not committed
    external/                    # not committed
    interim/                     # not committed
    processed/                   # not committed
    metadata/                    # small metadata TSVs can be committed if allowed

  src/
    ora/
      __init__.py
      io.py
      metadata.py
      qc.py
      aggregate.py
      gene_sets.py
      feature_engineering.py
      compositional.py
      trajectory.py
      age_model.py
      validation.py
      plotting.py
      reporting.py
      utils.py

  scripts/
    download_cellxgene.py
    inspect_h5ad.py
    build_sample_manifest.py
    aggregate_cell_counts.py
    aggregate_pseudobulk.py
    score_gene_sets.py
    build_feature_matrix.py
    run_compositional_models.py
    run_milo.R
    run_trajectory.py
    run_age_models.py
    run_external_validation.py
    make_figures.py
    make_report.py

  workflows/
    Snakefile
    rules/
      download.smk
      aggregate.smk
      features.smk
      models.smk
      validation.smk
      figures.smk

  notebooks/
    00_data_inspection.ipynb
    01_metadata_qc.ipynb
    02_cell_state_abundance.ipynb
    03_trajectory_pilot.ipynb
    04_ora_model_exploration.ipynb
    05_external_validation.ipynb
    99_manuscript_figures.ipynb

  results/
    tables/
    figures/
    models/
    reports/
    logs/

  tests/
    test_metadata.py
    test_aggregate.py
    test_features.py
    test_age_model.py
```

---

# Pipeline overview

## Stage 0: Data download and inspection

### Objective

Acquire Gateway data and inspect metadata without loading the full matrix into memory.

### Inputs

- CELLxGENE collection URL or dataset ID.
- Local output directory.

### Outputs

- `data/raw/gateway.h5ad`
- `results/reports/h5ad_schema.json`
- `data/metadata/gateway_obs_columns.tsv`
- `data/metadata/gateway_var_columns.tsv`
- `data/metadata/gateway_manifest.tsv`

### Implementation notes

Use `cellxgene_census.download_source_h5ad()` if dataset ID is available. If only the collection ID is available, first query the collection datasets or manually copy the dataset download URL from CELLxGENE. Keep a function that accepts either a dataset ID or direct H5AD URL.

Use AnnData backed mode:

```python
import anndata as ad
adata = ad.read_h5ad(path, backed="r")
print(adata.shape)
print(adata.obs.columns)
print(adata.var.columns)
print(adata.obsm.keys())
```

Do not call `adata.X.toarray()`.

---

## Stage 1: Cohort definition

### Objective

Define reproducible donor cohorts for discovery, sensitivity analyses, and exploratory disease anchoring.

### Primary discovery cohort

Healthy controls only.

Recommended filters:

- disease condition == healthy / control / non-NDD
- age not missing
- donor has at least 50 total olfactory neuronal lineage cells, or define separate thresholds by analysis
- donor has valid sex metadata
- donor has valid chemistry and collection method metadata

### Sensitivity cohorts

1. **All healthy controls**
2. **Healthy FLEX v2 only**
3. **Healthy device-only**
4. **Healthy donors with >= 50 mOSNs**
5. **Healthy donors with >= 500 neuronal lineage cells**
6. **Age-balanced subset across chemistry and collection method**, if possible

### Exploratory cohorts

- AD donors
- PD donors
- Combined NDD donors

These are never used to train ORA. They are projected after model freezing.

### Outputs

- `data/processed/cohort_manifest.tsv`
- `results/tables/cohort_summary.tsv`
- `results/figures/cohort_age_distribution.pdf`

---

## Stage 2: Donor-level cell-state composition

### Objective

Compute robust donor-level composition features.

### Inputs

- `adata.obs`
- cohort manifest

### Features

For each donor and sample:

1. Total cell count.
2. Total UMI count if available.
3. Cell counts by coarse cell type.
4. Cell counts by fine cell type.
5. Cell proportions by coarse and fine cell type.
6. Centered log-ratio transformed proportions.
7. Additive log-ratio transformed biologically motivated ratios.
8. Diversity metrics such as Shannon entropy across cell states.

### Key ratios

Use small pseudocounts and report sensitivity to pseudocount size.

```text
neuronal_fraction = (INP + iOSN + mOSN) / total_cells
mature_neuron_fraction = mOSN / total_cells
immature_to_mature_ratio = iOSN / (mOSN + epsilon)
progenitor_to_neuron_ratio = INP / (iOSN + mOSN + epsilon)
activated_to_quiescent_HBC = activated_HBC / (quiescent_HBC + epsilon)
stressed_to_fully_mature_mOSN = stressed_mOSN / (fully_mature_mOSN + epsilon)
respiratory_to_olfactory_epithelium = respiratory_epithelium / (olfactory_lineage + sustentacular + HBC + epsilon)
immune_to_epithelial_ratio = immune_cells / epithelial_cells
```

### Statistical models

Fit both continuous age and age-bin models.

Age as continuous:

```text
feature ~ ns(age, df=3) + sex + race_ethnicity + chemistry + collection_method + site + log10(total_cells)
```

If site is not available, use available technical covariates and run sensitivity by chemistry/method.

Age bins for visualization only:

```text
young: < 50
middle: 50-69
old: >= 70
```

Use continuous models for inference.

### Outputs

- `data/processed/donor_cell_state_counts.tsv`
- `data/processed/donor_cell_state_features.tsv`
- `results/tables/age_cell_state_associations.tsv`
- `results/figures/cell_state_by_age.pdf`
- `results/figures/lineage_ratio_by_age.pdf`

---

## Stage 3: Gene module scoring

### Objective

Compute biologically interpretable pathway and gene-program activity at cell-state and donor level.

### Strategy

Prefer cell-state-specific donor pseudobulk or cell-state mean expression features on a laptop. Full per-cell UCell scoring is useful, but not required for MVP if memory is limiting.

### Gene-set categories

#### Neurogenic lineage

- HBC identity: TP63, KRT5, KRT14, KRT15
- HBC activation / injury response: TP63, SOX2, EGFR, KRT6A, KRT6B, KRT17
- Progenitor / neuroblast: ASCL1, NEUROG1, NEUROD1, SOX2, MKI67, TOP2A
- Immature neuron: GAP43, STMN2, TUBB3, DCX
- Mature olfactory neuron: OMP, ADCY3, GNAL, GNG13, CNGA2, RTP1, RTP2
- Ciliogenesis / olfactory transduction: FOXJ1, IFT genes, BBS genes, ADCY3, CNGA2

These marker lists are starting points; verify and refine from dataset-specific markers.

#### Aging and stress

- Inflammatory response
- Interferon alpha/gamma response
- NF-kB signaling
- TNF signaling
- IL1 / IL6 / cytokine response
- Complement
- Cellular senescence
- DNA damage response
- Oxidative stress
- Proteostasis / unfolded protein response
- Ubiquitin-proteasome system
- Autophagy
- Lysosome / endosome
- Mitochondrial function / mitophagy

#### Neuronal maturity and function

- Axon development
- Synaptic signaling
- Dendrite / axon compartment
- Neuronal projection
- Ion channel activity
- Neurotransmission

### Methods

Recommended MVP:

1. Pseudobulk each donor x fine cell type.
2. Normalize pseudobulk counts with CPM/logCPM or DESeq2 variance-stabilized counts.
3. Score gene sets using mean z-score, ssGSEA, decoupler, or UCell on pseudobulk/cells.
4. Aggregate module scores to donor-level features.

Full version:

1. Run UCell per cell in chunks.
2. Average UCell scores by donor and fine cell type.
3. Fit age association models per cell type and module.

### Outputs

- `data/processed/donor_celltype_module_scores.tsv`
- `data/processed/donor_module_features.tsv`
- `results/tables/age_module_associations.tsv`
- `results/figures/module_heatmap_by_age.pdf`

---

## Stage 4: Pseudobulk differential expression with age

### Objective

Identify age-associated genes within major regenerative cell states.

### Cell states for DE

Primary:

- Quiescent HBC
- Activated HBC
- Early INP
- Late INP
- Early iOSN
- Late iOSN
- Early mature mOSN
- Fully mature mOSN
- Stressed mOSN
- Sustentacular cells
- Immune subsets if enough donors/cells

### Pseudobulk criteria

For a donor x cell-state aggregate:

- Minimum cells: start with >= 50 cells; sensitivity at >= 100 and >= 200.
- Minimum UMIs: use dataset-specific cutoffs.
- Minimum donors per model: at least 20 donors with valid aggregates.

### Model

Use edgeR/limma-voom, DESeq2, or dreamlet.

Continuous age model:

```text
counts ~ ns(age, df=3) + sex + race_ethnicity + chemistry + collection_method + log_lib_size
```

Alternative linear model:

```text
counts ~ age_scaled + sex + race_ethnicity + chemistry + collection_method + log_lib_size
```

Use linear age for interpretable effect sizes; use spline age for nonlinear discovery.

### Outputs

- `data/processed/pseudobulk_counts.h5ad` or matrices by cell type
- `results/tables/de_age_<celltype>.tsv`
- `results/tables/gsea_age_<celltype>.tsv`
- `results/figures/volcano_age_<celltype>.pdf`
- `results/figures/enrichment_dotplot_age.pdf`

---

## Stage 5: Differential abundance and neighborhood analysis

### Objective

Detect age-associated continuous cell-state shifts that cluster labels may miss.

### MVP method

Start with compositional models on annotated fine cell states.

Recommended methods:

- Dirichlet-multinomial regression.
- scCODA.
- `propeller` from speckle/limma ecosystem.
- Beta-binomial regression for individual lineages.

### Full method

Run Milo on lineage subsets.

Subsets:

1. Basal + neuronal lineage cells: HBC, INP, iOSN, mOSN.
2. Epithelial support cells: sustentacular, microvillous, Bowman gland, respiratory epithelium.
3. Immune cells.

Model:

```text
nhood_abundance ~ age_scaled + sex + race_ethnicity + chemistry + collection_method + log10(total_cells)
```

If enough donors per category:

```text
nhood_abundance ~ ns(age, df=3) + sex + chemistry + collection_method + log10(total_cells)
```

### Computational strategy for Milo

Laptop mode:

- Use existing embedding if present, e.g. X_scANVI.
- Restrict to lineage cells.
- Downsample to max 1,000 to 3,000 cells per donor per broad lineage.
- Build KNN on 30 to 50 dimensions.

External compute mode:

- Use all lineage cells passing QC.
- Run separate Milo objects for each lineage compartment.

### Outputs

- `results/tables/milo_age_lineage.tsv`
- `results/tables/milo_age_immune.tsv`
- `results/figures/milo_da_umap_age.pdf`
- `results/figures/milo_nhood_effects_by_pseudotime.pdf`

---

## Stage 6: Trajectory and bottleneck analysis

### Objective

Create an interpretable cross-sectional map of neurogenic progression and identify age-associated bottleneck regions.

### Cells included

- HBC subtypes
- INP subtypes
- iOSN subtypes
- mOSN subtypes

Exclude immune and respiratory cells from the primary trajectory.

### Embedding

Preferred:

- Use Gateway-provided scANVI latent embedding if included.

Alternative:

- Compute PCA on highly variable genes using a balanced cell subset.
- Use Harmony or scVI only if needed.

### Pseudotime methods

MVP:

- Diffusion pseudotime or PAGA in Scanpy.
- Root at quiescent HBC / activated HBC.
- Terminal states: fully mature mOSN and stressed mOSN.

Optional:

- Palantir for branch probabilities.
- Optimal transport-inspired density movement only if multiple age bins are treated as pseudo-time cohorts.
- Avoid RNA velocity unless spliced/unspliced layers are confirmed available and reliable for 10x FLEX.

### Bottleneck metrics

Define pseudotime bins, e.g. 100 bins across lineage pseudotime. For each donor:

```text
density[d, bin] = number of donor cells in bin / total donor lineage cells
```

Fit:

```text
density_bin ~ age_scaled + sex + chemistry + collection_method + log10(total_lineage_cells)
```

Identify bins where age effect is most negative or positive after FDR correction.

Transition ratios:

```text
HBC_activation_ratio = activated_HBC / (quiescent_HBC + epsilon)
HBC_to_INP_ratio = INP / (activated_HBC + epsilon)
INP_to_iOSN_ratio = iOSN / (INP + epsilon)
iOSN_to_mOSN_ratio = fully_mature_mOSN / (iOSN + epsilon)
stressed_mOSN_ratio = stressed_mOSN / (fully_mature_mOSN + epsilon)
```

Interpretation:

- These are cross-sectional density proxies, not direct measured cell flux.
- Use conservative language: "candidate bottleneck," "lineage-density shift," "maturation imbalance."

### Outputs

- `data/processed/lineage_pseudotime.tsv`
- `data/processed/donor_pseudotime_density.tsv`
- `results/tables/pseudotime_age_associations.tsv`
- `results/figures/lineage_pseudotime_umap.pdf`
- `results/figures/pseudotime_density_by_age.pdf`
- `results/figures/bottleneck_summary.pdf`

---

## Stage 7: ORA model

### Objective

Train a donor-level model that predicts chronological age from interpretable single-cell features and uses the model to define olfactory regenerative age.

### Input feature matrix

Rows: donors.

Columns:

1. Cell-state CLR proportions.
2. Biological lineage ratios.
3. Cell-state-specific module scores.
4. Pseudotime density features.
5. Pseudobulk age-associated gene-program scores.
6. Optional technical features only for adjustment, not biological interpretation.

### Exclusions

Do not include raw total cell count, sequencing depth, chemistry, collection method, or site as predictive biological features in the main ORA model. Use them as covariates for residualization or sensitivity analysis.

### Feature preprocessing

1. Remove features missing in > 30% of donors.
2. Impute remaining missing values with training-fold medians.
3. Standardize features within training fold.
4. Residualize biological features against technical covariates within training fold, or include technical covariates only in a post-hoc calibration model.
5. Use donor-level train/test splits only.

### Models

Run at least four models:

1. **Null model:** mean age or covariates-only model.
2. **Elastic net:** primary interpretable model.
3. **Gradient boosting:** XGBoost or LightGBM for nonlinear effects.
4. **Random forest:** robustness baseline.

Optional:

- Generalized additive model for smooth effects.
- Bayesian regression for uncertainty.
- Stacked model if justified, but keep final model simple for publication.

### Cross-validation

Use nested donor-level cross-validation.

Outer loop:

- 5 folds, stratified by age bin, sex, chemistry, and collection method if possible.

Inner loop:

- Hyperparameter tuning.

Repeated CV:

- 10 repeats if dataset size permits.

Never split cells independently; that leaks donor identity.

### Metrics

- Mean absolute error (MAE)
- Root mean squared error (RMSE)
- R-squared
- Spearman correlation between predicted and chronological age
- Calibration slope
- Permutation p-value by shuffling age labels within major technical strata

### ORA definition

```text
ORA_donor = out-of-fold predicted chronological age from biological features
```

### ORAA definition

Fit on healthy donors only:

```text
ORA_donor ~ chronological_age + sex + race_ethnicity + chemistry + collection_method + site
```

Then:

```text
ORAA = residual from this model
```

Alternative:

```text
ORAA = ORA - expected_ORA_given_age_and_covariates
```

### NCI definition

Use the sign of model features to define a biological neurogenic capacity score.

Example:

```text
NCI = z(beneficial regenerative features) - z(aging/stress features)
```

Where beneficial features may include progenitor abundance, immature neuron abundance, healthy mature neuron abundance, synaptic/ciliary maturation scores, and low inflammatory stress. The final included features should be derived from the trained elastic-net coefficients and biological interpretation.

### Interpretability

Use:

- Elastic-net coefficients.
- SHAP values for XGBoost.
- Feature stability across folds.
- Partial dependence plots for top features.
- Grouped feature importance by biological category.

### Outputs

- `data/processed/ora_feature_matrix.tsv`
- `results/models/ora_elasticnet.joblib`
- `results/models/ora_xgboost.joblib`
- `results/tables/ora_model_performance.tsv`
- `results/tables/ora_feature_importance.tsv`
- `results/tables/donor_ora_scores.tsv`
- `results/figures/ora_predicted_vs_chronological_age.pdf`
- `results/figures/ora_feature_importance.pdf`
- `results/figures/ora_shap_summary.pdf`

---

## Stage 8: External validation

### Objective

Show that ORA captures biologically meaningful olfactory aging and disease-relevant biology beyond the Gateway cohort.

### Validation mode A: gene-set transfer

Use the top positive and negative ORA genes/programs as gene sets.

For each external dataset:

1. Map genes to common symbols or Ensembl IDs.
2. Score ORA-positive and ORA-negative modules in relevant cell types.
3. Compare by published condition:
   - Presbyosmia vs normosmia.
   - AD stage vs control.
   - Preclinical AD vs control.
4. Use donor-level aggregation where donor metadata are available.

### Validation mode B: feature transfer

If external datasets have compatible cell labels:

1. Harmonize cell types to Gateway labels.
2. Compute cell-state proportions and module scores.
3. Apply frozen ORA model if feature overlap is sufficient.
4. Report lower-confidence because cell type labels and chemistry differ.

### Validation mode C: overlap enrichment

If raw external data are difficult to process quickly:

1. Extract published DE gene lists or supplementary tables.
2. Test overlap with ORA-associated genes using Fisher's exact test and rank-based GSEA.
3. Report directionality where possible.

### External validation hypotheses

1. Presbyosmic HBCs show increased ORA-positive inflammatory/stress modules.
2. Presbyosmic HBCs show altered TP63/differentiation-blockade modules.
3. AD olfactory biopsies show increased ORA-positive immune/neuronal injury modules.
4. ORA-positive genes overlap with AD olfactory neuron inflammatory injury programs more than random age-associated genes.

### Outputs

- `results/tables/external_validation_summary.tsv`
- `results/tables/ora_overlap_oliva.tsv`
- `results/tables/ora_overlap_danniballe.tsv`
- `results/figures/external_validation_heatmap.pdf`
- `results/figures/ora_modules_external_boxplots.pdf`

---

## Stage 9: Exploratory NDD projection

### Objective

Use Gateway AD and PD donors only after the ORA model is frozen.

### Analysis

1. Compute ORA and ORAA for AD and PD donors.
2. Compare NDD donors to age-matched healthy donors.
3. Use nearest-neighbor matching on age, sex, chemistry, and collection method if possible.
4. Report effect sizes with wide confidence intervals.
5. Treat results as hypothesis-generating.

### Model

```text
ORAA ~ disease_status + age + sex + chemistry + collection_method
```

Use caution due to small NDD n.

### Outputs

- `results/tables/ndd_ora_projection.tsv`
- `results/figures/ndd_projection_ora.pdf`

---

# Manuscript figure plan

## Figure 1: Study design and regenerative aging framework

Panels:

A. Schematic: olfactory neurogenic lineage and ORA concept.
B. Gateway healthy cohort age distribution.
C. Coarse/fine cell-state composition across donors.
D. Analysis pipeline overview.

Main message:

The Gateway atlas supports a donor-level model of olfactory regenerative aging.

## Figure 2: Age-associated cell-state composition

Panels:

A. Cell-state proportions across age.
B. Differential abundance effect-size heatmap.
C. Key lineage ratios by age.
D. Sensitivity across all healthy, FLEX v2 only, device-only.

Main message:

Aging is associated with coordinated shifts in regenerative, neuronal, epithelial, and immune compartments.

## Figure 3: Neurogenic trajectory and bottleneck

Panels:

A. UMAP or diffusion map of HBC -> INP -> iOSN -> mOSN lineage.
B. Pseudotime marker validation.
C. Donor density along pseudotime by age.
D. Age-effect curve across pseudotime bins.
E. Bottleneck transition ratio summary.

Main message:

Aging localizes to specific lineage-density transitions, suggesting a candidate regenerative bottleneck.

## Figure 4: ORA model

Panels:

A. Out-of-fold predicted vs chronological age.
B. Model performance vs null.
C. Top interpretable features.
D. Feature groups contributing to ORA.
E. ORAA distribution.

Main message:

Donor age is predictable from interpretable regenerative single-cell features.

## Figure 5: Gene programs underlying ORA

Panels:

A. Module-age association heatmap by cell type.
B. HBC/progenitor inflammatory and differentiation programs.
C. Neuronal maturation/synaptic/ciliary/proteostasis programs.
D. GSEA dot plot for age-associated pseudobulk DE.

Main message:

ORA links aging to inflammation, stem-cell stress, and reduced neuronal maturation/homeostasis programs.

## Figure 6: External validation and disease anchoring

Panels:

A. ORA-positive modules in Oliva presbyosmia dataset.
B. ORA-positive modules in D'Anniballe AD olfactory biopsy dataset.
C. Overlap between ORA genes and published presbyosmia/AD signatures.
D. Exploratory Gateway AD/PD ORAA projection.

Main message:

ORA-associated programs generalize to independent olfactory aging/disease datasets.

---

# Statistical analysis details

## Covariates

Always inspect correlation between age and:

- chemistry
- collection method
- site
- sex
- race / ethnicity
- total cells
- neuronal yield
- sequencing depth
- disease status

### Main covariate set

```text
age + sex + race_ethnicity + chemistry + collection_method + log10(total_cells)
```

### Extended covariate set

```text
age + sex + race_ethnicity + chemistry + collection_method + site + log10(total_cells) + log10(total_UMI)
```

Use extended model only if site and total UMI are available and not too collinear.

## Multiple testing

- Use Benjamini-Hochberg FDR for cell-state, module, gene, and neighborhood tests.
- Report exact number of tests per family.
- Keep primary hypotheses separate from exploratory analyses.

## Effect sizes

Report effect sizes alongside p-values:

- Change per 10 years of age.
- Old vs young contrast estimated from continuous model.
- Standardized beta.
- Cohen's d for visualization only.

## Sensitivity analyses

Minimum required:

1. Healthy-only discovery.
2. FLEX v2-only healthy controls.
3. Device-only healthy controls.
4. Downsample equal cells per donor.
5. Exclude donors with very low olfactory neuron yield.
6. Include/exclude race/ethnicity depending on missingness.
7. Permute age labels within chemistry/method strata.
8. Refit ORA without direct cell-yield features.

---

# Minimum viable paper (MVP)

The fastest publishable version can be completed without full-cell deep learning.

## MVP analyses

1. Download H5AD and extract metadata.
2. Healthy donor cohort definition.
3. Donor-level cell-state composition and lineage ratios.
4. Pseudobulk module scores by donor x cell type.
5. Age association models for cell states and modules.
6. Elastic-net ORA model with nested CV.
7. Downsampled trajectory with pseudotime density.
8. External validation using Oliva and/or D'Anniballe published data or supplementary gene lists.
9. Manuscript figures 1-6.

## MVP novelty

The novelty is the **donor-level regenerative aging model** and **bottleneck framing**, not the use of a new neural network.

## MVP risk

If ORA predicts age weakly, the paper can still work if:

- Age-associated cell-state changes are robust.
- The bottleneck analysis is clear.
- External validation supports presbyosmia/AD relevance.

But if both age effects and validation are weak, pivot to a methods note / resource reanalysis.

---

# Stretch goals

## Stretch 1: Full Milo neighborhood analysis

Impact: high.
Compute: external recommended.

## Stretch 2: cNMF program discovery

Impact: medium-high.
Compute: external recommended for large subsets.

## Stretch 3: Transfer ORA to external datasets directly

Impact: high if it works.
Risk: batch and label differences.

## Stretch 4: Integrate SEA-AD brain aging modules

Impact: high.
Risk: scope creep.

## Stretch 5: Build a public web explorer

Impact: medium.
Risk: distracts from paper.

---

# Engineering milestones

## Week 1: Data ingestion and cohort setup

Deliverables:

- H5AD downloaded.
- Metadata schema report.
- Healthy cohort manifest.
- Cell-state count table.
- Initial cohort QC figures.

Codex tasks:

1. Create repo structure.
2. Implement `download_cellxgene.py`.
3. Implement `inspect_h5ad.py`.
4. Implement `build_sample_manifest.py`.
5. Implement `aggregate_cell_counts.py`.

## Week 2: Feature engineering and age associations

Deliverables:

- Donor-level feature matrix.
- Cell-state age models.
- Module score pipeline.
- First results figures.

Codex tasks:

1. Implement compositional transforms.
2. Implement biological ratio features.
3. Implement gene-set loading.
4. Implement donor x cell-type module scoring.
5. Implement statsmodels age association scripts.

## Week 3: ORA model

Deliverables:

- Nested CV ORA model.
- Performance tables.
- Feature importance.
- ORAA scores.

Codex tasks:

1. Implement sklearn pipelines.
2. Implement nested CV splitter stratified by age/sex/chemistry/method.
3. Implement null permutation tests.
4. Implement SHAP or coefficient stability.
5. Implement model report.

## Week 4: Trajectory and bottleneck analysis

Deliverables:

- Downsampled lineage trajectory.
- Pseudotime density by donor.
- Bottleneck association table.
- Figure 3.

Codex tasks:

1. Implement lineage subsetter.
2. Implement balanced downsampling.
3. Implement pseudotime wrapper.
4. Implement donor-density features.
5. Implement bottleneck plots.

## Week 5: External validation

Deliverables:

- External validation data loaded or gene lists curated.
- ORA module transfer results.
- Figure 6.

Codex tasks:

1. Implement external dataset import adapters.
2. Implement gene ID harmonization.
3. Implement module scoring on external datasets.
4. Implement overlap/GSEA tests.
5. Implement validation report.

## Week 6: Manuscript package

Deliverables:

- Final figures.
- Final tables.
- Reproducible pipeline run.
- Draft methods text.
- Draft results outline.

Codex tasks:

1. Implement `make_figures.py`.
2. Implement `make_report.py`.
3. Write methods from config and outputs.
4. Freeze environment.
5. Produce reproducibility checklist.

---

# Codex implementation tickets

## Ticket 1: Repository scaffold

Acceptance criteria:

- Repo follows structure above.
- `make setup` creates environment.
- `make test` runs unit tests.
- README has quickstart.

## Ticket 2: H5AD inspector

Acceptance criteria:

- Works in backed mode.
- Exports obs/var/obsm schema.
- Prints memory-safe summary.
- Does not load full matrix.

## Ticket 3: Cohort manifest builder

Acceptance criteria:

- User can map actual metadata column names in YAML.
- Outputs donor/sample manifest.
- Flags missing age, disease, sex, chemistry, collection method.
- Generates cohort summary table.

## Ticket 4: Cell-state aggregation

Acceptance criteria:

- Counts cells by donor, sample, coarse cell type, fine cell type.
- Computes proportions and CLR transforms.
- Handles zero counts with pseudocount.
- Unit tests on toy data.

## Ticket 5: Gene-set module scorer

Acceptance criteria:

- Loads GMT/CSV gene sets.
- Maps symbols to Ensembl or vice versa.
- Scores pseudobulk and/or cell chunks.
- Aggregates to donor x cell-type table.

## Ticket 6: ORA feature matrix builder

Acceptance criteria:

- Merges composition, ratios, module scores, pseudotime density when available.
- Applies missingness filters.
- Writes feature dictionary with feature categories.

## Ticket 7: ORA model trainer

Acceptance criteria:

- Donor-level nested CV.
- Elastic net, random forest, XGBoost/LightGBM, null model.
- Outputs out-of-fold predictions.
- Outputs performance metrics and feature importances.
- Implements permutation testing.

## Ticket 8: Trajectory pipeline

Acceptance criteria:

- Selects lineage cells.
- Balanced downsampling by donor.
- Computes pseudotime.
- Exports donor pseudotime density.
- Generates lineage marker plots.

## Ticket 9: Differential abundance

Acceptance criteria:

- Runs cluster-level compositional models.
- Optionally runs Milo through R script.
- Outputs FDR-adjusted tables.

## Ticket 10: External validation

Acceptance criteria:

- Imports at least one external dataset or published gene list.
- Harmonizes genes.
- Scores ORA modules.
- Outputs validation statistics and plots.

## Ticket 11: Figure generation

Acceptance criteria:

- Recreates all manuscript figures from processed outputs.
- Saves PDF and PNG.
- Uses consistent labels and legends.

## Ticket 12: Report generation

Acceptance criteria:

- Generates HTML or Markdown report with all key outputs.
- Includes package versions and config.
- Includes cohort sizes and exclusions.

---

# Example config files

## `configs/gateway.yaml`

```yaml
collection_id: "8b35aa1f-6bcf-4a51-abc3-a3f336a44ae6"
dataset_id: null
h5ad_path: "data/raw/gateway.h5ad"

columns:
  donor_id: "donor_id"
  sample_id: "sample_id"
  age: "age"
  sex: "sex"
  race_ethnicity: "race_ethnicity"
  disease: "disease_condition"
  chemistry: "flex_chemistry_version"
  collection_method: "collection_method"
  coarse_cell_type: "coarse_cell_type"
  fine_cell_type: "fine_cell_type"
  n_counts: "nCount_RNA"
  n_genes: "nFeature_RNA"
  percent_mito: "percent_mito"

healthy_values:
  - "healthy"
  - "control"
  - "cognitively normal"

lineage_cell_types:
  - "Quiescent HBC"
  - "Activated HBC"
  - "Cycling HBC"
  - "Early INP"
  - "Late INP"
  - "Early iOSN"
  - "Late iOSN"
  - "Early mature mOSN"
  - "Fully mature mOSN"
  - "Stressed mOSN"
```

## `configs/models.yaml`

```yaml
random_seed: 42
outer_cv_folds: 5
outer_cv_repeats: 10
inner_cv_folds: 5
age_bins:
  young: [0, 49]
  middle: [50, 69]
  old: [70, 120]

missingness_max_fraction: 0.30
imputation: "median"
scaling: "standard"

models:
  elastic_net:
    alphas: [0.001, 0.01, 0.1, 1.0, 10.0]
    l1_ratios: [0.1, 0.5, 0.9, 1.0]
  random_forest:
    n_estimators: [500]
    max_depth: [3, 5, 10, null]
  xgboost:
    n_estimators: [200, 500]
    max_depth: [2, 3, 5]
    learning_rate: [0.01, 0.05, 0.1]
```

---

# Example command-line interface

```bash
make setup
make download-gateway
make inspect
make cohort
make aggregate
make features
make model-ora
make trajectory
make validate
make figures
make report
```

Equivalent direct calls:

```bash
python scripts/inspect_h5ad.py --config configs/gateway.yaml
python scripts/build_sample_manifest.py --config configs/gateway.yaml
python scripts/aggregate_cell_counts.py --config configs/gateway.yaml
python scripts/build_feature_matrix.py --gateway configs/gateway.yaml --features configs/features.yaml
python scripts/run_age_models.py --features data/processed/ora_feature_matrix.tsv --config configs/models.yaml
```

---

# Key tables to produce

## Table 1: Cohort summary

Columns:

- cohort
- donors
- samples
- cells
- median age
- age IQR
- sex counts
- chemistry counts
- collection method counts
- major cell-type counts

## Table 2: Age-associated cell states

Columns:

- cell_state
- broad_lineage
- beta_per_10_years
- standard_error
- p_value
- fdr
- direction
- sensitivity_all_healthy
- sensitivity_flex_v2
- sensitivity_device_only

## Table 3: ORA model performance

Columns:

- model
- MAE
- RMSE
- R2
- Spearman_r
- calibration_slope
- permutation_p

## Table 4: ORA feature importance

Columns:

- feature
- feature_category
- cell_type
- coefficient_or_importance
- stability_across_folds
- direction
- biological_interpretation

## Table 5: ORA-associated pathways

Columns:

- cell_type
- pathway
- beta_per_10_years
- p_value
- fdr
- direction
- interpretation

## Table 6: External validation

Columns:

- dataset
- comparison
- cell_type
- ORA_module
- effect_size
- p_value
- fdr
- direction_consistent

---

# Risks and mitigation

## Risk 1: Age is confounded with chemistry or collection method

Mitigation:

- Report age distribution by chemistry and method.
- Use covariate adjustment.
- Run FLEX v2-only and device-only sensitivity analyses.
- Use permutation within technical strata.
- Do not interpret features that disappear in all technical sensitivity analyses.

## Risk 2: Cell yield drives apparent aging signal

Mitigation:

- Include total cells and lineage cells as covariates.
- Remove raw yield features from ORA model.
- Downsample equal cells per donor.
- Use proportions and CLR features.

## Risk 3: Disease donors are older and all device/FLEX v2

Mitigation:

- Train ORA only on healthy donors.
- Use NDD only as exploratory projection.
- Use age-matched controls and wide confidence intervals.

## Risk 4: External datasets differ in chemistry and cell labels

Mitigation:

- Validate gene programs, not exact ORA values, unless feature transfer is robust.
- Use broad cell types for external validation.
- Report gene overlap and module scores separately.

## Risk 5: Trajectory analysis overclaims lineage flux

Mitigation:

- Use language like "cross-sectional lineage density" and "candidate bottleneck."
- Support with marker gradients and independent abundance ratios.
- Avoid claims about actual rates of differentiation without longitudinal data.

---

# Draft methods text skeleton

## Cohort definition

We downloaded the Gateway human olfactory epithelium atlas from CELLxGENE and analyzed healthy control donors with non-missing chronological age as the primary discovery cohort. Donors with neurodegenerative disease were excluded from model training and reserved for exploratory projection. Metadata fields were harmonized to donor ID, sample ID, age, sex, race/ethnicity, disease status, FLEX chemistry, collection method, and coarse/fine cell-state labels. All analyses were performed at donor level unless otherwise noted.

## Cell-state abundance analysis

For each donor, we computed cell counts and proportions across Gateway fine cell states. To account for the compositional nature of cell proportions, proportions were transformed using centered log-ratio transformation with a fixed pseudocount. We fit linear models testing association between each transformed cell-state abundance and chronological age, adjusting for sex, race/ethnicity, chemistry, collection method, and total cell count. Multiple testing was controlled using Benjamini-Hochberg FDR.

## Gene-program scoring

Gene sets representing olfactory neurogenesis, neuronal maturation, inflammation, proteostasis, lysosomal biology, mitochondrial stress, senescence, cilia, and synaptic function were curated from MSigDB, Gene Ontology, and literature marker sets. Scores were computed for donor x cell-state pseudobulks and aggregated into donor-level features. Age associations were tested using covariate-adjusted linear models.

## ORA model

We trained donor-level models to predict chronological age from biological single-cell features, including cell-state composition, lineage ratios, cell-state-specific module scores, and pseudotime-density features. Model performance was evaluated using nested donor-level cross-validation. The primary ORA model used elastic-net regression for interpretability, with random forest and gradient-boosted tree models as nonlinear benchmarks. Olfactory Regenerative Age was defined as the out-of-fold predicted age. Olfactory Regenerative Age Acceleration was defined as the residual ORA after adjustment for chronological age and technical covariates.

## Trajectory analysis

Cells from the basal-to-neuronal lineage were subset and embedded using the available latent representation or PCA on highly variable genes. Pseudotime was inferred with root states defined by HBC markers and terminal states defined by mature olfactory sensory neuron markers. Donor-level density across pseudotime bins was computed and tested for age association using covariate-adjusted models. Regions with significant age-associated density loss or gain were interpreted as candidate regenerative bottlenecks.

## External validation

ORA-associated gene programs were evaluated in independent human olfactory datasets. For each dataset, genes were harmonized to common identifiers, module scores were computed in comparable cell populations, and differences across published clinical groups were tested at donor level where metadata were available. Published gene signatures were also tested for overlap with ORA-associated genes using rank-based enrichment and Fisher exact tests.

---

# Draft results claims to aim for

Use these only if supported by the actual analysis.

1. Healthy human olfactory epithelium shows coordinated age-associated remodeling across basal, neuronal, epithelial, and immune compartments.
2. Age-associated effects are not uniform across the neurogenic lineage but concentrate at a specific maturation transition.
3. Donor chronological age can be predicted from interpretable regenerative features, defining ORA.
4. ORA-positive features include inflammatory/stress programs and stressed mature neuronal states.
5. ORA-negative features include neurogenic progenitor activity, neuronal maturation, cilia/olfactory transduction, synaptic signaling, or proteostasis programs.
6. ORA-associated programs overlap with independent presbyosmia and AD olfactory biopsy signatures.
7. Gateway AD/PD donors show exploratory shifts toward older ORA profiles, but this requires larger validation.

---

# Journal targets

## Fast and realistic

- Aging Cell
- GeroScience
- eLife
- iScience
- Cell Reports
- Communications Biology
- npj Aging
- npj Systems Biology and Applications

## Higher ambition

- Cell Genomics
- Genome Biology
- Nature Aging
- Nature Communications
- Genome Medicine

## ML/methods framing

- Bioinformatics
- PLOS Computational Biology
- Genome Biology
- Nature Computational Science, only if the method becomes substantially novel beyond standard models

---

# Decision rules for pivoting

## Continue with full ORA paper if:

- ORA model beats null with convincing donor-level CV.
- At least several cell states or modules associate with age after covariate adjustment.
- Sensitivity analyses do not erase the main pattern.
- External validation supports at least one core biological module.

## Pivot to biology-only aging paper if:

- ORA prediction is weak, but abundance/trajectory/module age effects are strong.

New title:

**Population-scale single-cell analysis identifies regenerative bottlenecks in aging human olfactory neuroepithelium**

## Pivot to methods/resource note if:

- Main effects are modest but pipeline and benchmarking are strong.

New title:

**A reproducible framework for donor-level aging analysis in large single-cell atlases of accessible human neural tissue**

---

# Practical first 48 hours

1. Download or start download of the Gateway H5AD.
2. Inspect `.obs` and `.obsm` columns.
3. Confirm exact labels for age, disease, cell type, chemistry, and collection method.
4. Count donors and cells after filtering to healthy controls with non-missing age.
5. Make the first donor-level cell-state proportion table.
6. Plot major lineage proportions against age.
7. Run a quick elastic-net model on cell-state proportions only.
8. Decide whether signal is strong enough to proceed.

Fast sanity checks:

```text
Does mature neuron fraction change with age?
Does stressed mOSN fraction change with age?
Does activated/quiescent HBC ratio change with age?
Does immune/epithelial ratio change with age?
Can cell-state proportions predict age above null?
Are age effects preserved within FLEX v2 only?
```

If the answer to at least three of these is yes, the project is likely worth pushing to preprint.
