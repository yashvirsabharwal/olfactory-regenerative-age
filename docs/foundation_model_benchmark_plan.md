# Foundation-Model and Aging-Clock Benchmark Plan

Status: active benchmark plan  
Created: 2026-06-25  
Primary target: Gate B computational/SOTA submission

## Objective

Benchmark ORA against modern single-cell representation and aging-clock-style baselines under the same donor-level prediction problem used by the core ORA analyses. The goal is not to rebrand ORA as a universal biological clock. The goal is to show whether general pretrained transcriptomic representations, olfactory-specific interpretable features, or hybrids provide the best donor-level healthy-aging signal in human olfactory epithelium.

## Benchmark Principles

- Donor is the unit of prediction and splitting.
- Chronological age is predicted only for healthy, age-known donors eligible for ORA training.
- All model families must use the exported donor fold table unless a model-specific adaptation is explicitly logged.
- Cell-level embeddings are aggregated to donor-level features before age modeling.
- Performance is reported as MAE, RMSE, R2, Spearman rho, calibration slope, runtime, memory, gene coverage, and interpretability tier.
- Foundation-model failures are useful results if the failure reason is concrete: unavailable checkpoint, incompatible gene vocabulary, infeasible memory, unsupported local hardware, or license/access restrictions.

## Models

| Model family | Why include it | Planned representation | Current status |
| --- | --- | --- | --- |
| ORA composition | Native interpretable baseline | Donor-level proportions, CLR features, and lineage ratios | Complete |
| ORA modules | Biological program baseline | Donor-level curated module scores | Complete |
| scVI | Modern atlas latent baseline | Donor-level global and cell-state-specific latent summaries from full 4M scVI | Complete |
| ORA+scVI hybrid | Interpretable-plus-latent benchmark | Composition, modules, and donor scVI embeddings | Complete |
| Geneformer | Transformer pretrained across broad human single-cell transcriptomes | Rank-value tokenized cells, pooled cell embeddings, donor aggregation | Planned |
| scGPT | Generative single-cell foundation model with whole-human checkpoint | Cell embeddings from whole-human or continual-pretrained checkpoint, donor aggregation | Planned |
| scFoundation | Large single-cell transcriptomic foundation model | Cell embeddings through local model code or CLI/API if local weights are blocked | Planned |
| Aging-clock-style expression baseline | Reviewer-facing comparison to age-clock expectations | Donor-level expression PCs or fold-internal expression embeddings trained under ORA donor splits | Planned |

## External Model Context

Geneformer is a transformer model family pretrained on large human single-cell transcriptome corpora. The current Hugging Face model card describes V1 trained on about 30 million transcriptomes and V2 models trained on about 104 million transcriptomes, with rank-value transcriptome encoding and GPU recommended for efficient use: https://huggingface.co/ctheodoris/Geneformer

scGPT is the official codebase for "Towards Building a Foundation Model for Single-Cell Multi-omics Using Generative AI." Its repository recommends the whole-human checkpoint for most applications and lists pretrained models including whole-human, brain, blood, heart, lung, kidney, and pan-cancer variants. It also notes version-sensitive dependencies such as flash-attention: https://github.com/bowang-lab/scGPT

scFoundation is described by the project repository as a 100M-parameter model trained on over 50 million human single-cell transcriptomes, with model weights/code and downstream task examples. The old API was discontinued in 2024 and the repository points users to a newer platform/CLI path, so local-weight feasibility must be checked before treating it as a runnable baseline: https://github.com/biomap-research/scFoundation

## Input Subsets

The benchmark subset builder writes three fixed H5AD inputs from `data/raw/gateway.h5ad`:

- `data/processed/foundation_benchmark_lineage_subset.h5ad`: olfactory basal-to-neuronal lineage cells enriched for HBC, INP, iOSN, and mOSN states.
- `data/processed/foundation_benchmark_epithelial_subset.h5ad`: broad epithelial subset including olfactory, respiratory, sustentacular, secretory, glandular, and neuronal states.
- `data/processed/foundation_benchmark_allcell_subset.h5ad`: all-cell donor/fine-cell-type stratified benchmark subset.

Default cell caps are 120,000 lineage cells, 180,000 epithelial cells, and 250,000 all cells. These defaults are intentionally modest enough for repeated tokenizer passes on a 32 GB laptop, but can be raised on CPU/GPU infrastructure with `FOUNDATION_LINEAGE_CELLS`, `FOUNDATION_EPITHELIAL_CELLS`, and `FOUNDATION_ALLCELL_CELLS`.

The builder also writes:

- `results/tables/foundation_benchmark_subset_manifest.tsv`
- `results/tables/foundation_benchmark_donor_splits.tsv`
- `results/tables/foundation_benchmark_gene_manifest.tsv`
- `resources/foundation_benchmark/gateway_gene_symbols.txt`
- `resources/foundation_benchmark/gateway_gene_ids.txt`

## Donor-Level Aggregation

For each foundation model, cell embeddings will be summarized into donor-level features using the same pattern already implemented for scVI:

- global donor mean and standard deviation per embedding dimension;
- optional donor-by-fine-cell-type mean embedding dimensions for cell states with enough cells;
- QC table recording cells, donors, gene coverage, missing embedding fraction, runtime, and peak memory;
- identical repeated donor-level CV age models used by ORA and scVI baselines.

## Failure Criteria

A model is marked `blocked` or `deferred` rather than silently skipped if any of the following applies:

- checkpoint cannot be obtained under usable license/access conditions;
- required tokenizer vocabulary cannot map at least 50 percent of Gateway genes or a prespecified core marker/module gene panel;
- local CPU/MPS runtime cannot process the 120k lineage subset within a practical dry-run window;
- GPU-only dependency cannot be installed locally and no remote GPU run has been staged;
- model produces embeddings but donor-level aggregation has fewer than 40 eligible healthy training donors.

## Claim Language

If a foundation model outperforms ORA and ORA+scVI, claim that generic pretrained transcriptomic representations capture stronger donor-level age signal than the current interpretable ORA feature set, then reposition ORA as the interpretable biological decomposition.

If ORA+scVI outperforms foundation models, claim that an olfactory-aware interpretable-plus-latent representation is competitive with or stronger than generic single-cell foundation embeddings for this donor-level tissue-aging task.

If foundation models fail or perform similarly, claim that the analysis benchmarks against current representation-learning approaches where feasible and transparently documents model-access, vocabulary, and compute limitations.

No result unlocks an absolute biological-age-clock claim without stronger external validation and independent donor cohorts.
