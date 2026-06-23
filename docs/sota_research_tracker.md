# SOTA Research Tracker

Updated: 2026-06-23

Purpose: track the research work needed to make ORA a state-of-the-art olfactory regenerative aging resource, not merely a submission-ready manuscript. This is the project command center for scientific upgrades, validation, mechanism, and product-grade atlas deliverables.

Working thesis: the current project supports a modest, reproducible, interpretable healthy-donor olfactory regenerative aging axis. The SOTA target is stronger: a validated cellular, regulatory, and tissue-context model of how human olfactory epithelial regeneration changes with age.

## Status Key

| Status | Meaning |
| --- | --- |
| complete | Evidence exists, artifact is present or intentionally documented, and claim language is locked. |
| in progress | Work has started and the next artifact is clear. |
| needs work | Scientifically important and not yet started or not yet sufficient. |
| blocked | Requires unavailable data, compute, tissue, collaborator action, or external response. |
| stretch | High-upside work that would make the paper substantially stronger but may exceed the first submission. |

## Current Project State

| Area | Current state | SOTA gap | Decision |
| --- | --- | --- | --- |
| Core ORA model | Donor-split healthy model, repeated CV, shuffled-age null, calibration, sensitivity, and feature interpretation are present. | Age prediction is modest and under-dispersed, so it should not be sold as an age clock. | Keep as a reproducible regenerative tissue-state axis. |
| Atlas scale | Gateway provides a very large olfactory epithelial single-cell resource with full-scale scVI and Milo-style outputs. | Need stronger cross-cohort harmonization, foundation-model baselines, and cell-state transfer benchmarks. | Build an integrated external atlas and benchmark table. |
| Latent neighborhoods | Full 4M scVI/Milo-style summaries, edgeR parity, official MiloR subset sensitivity, and age-bin checks are present. | Reviewer will still ask whether neighborhoods represent fate, density, technical structure, or biology. | Treat as neighborhood remodeling unless fate/spatial/perturbation evidence is added. |
| External validation | GSE184117 and contextual sources are small or label-limited; final-search logging exists. | This is the largest scientific weakness. | Push a full validation atlas, request labels, and classify each dataset by validation strength. |
| Disease projection | AD/PD projection exists but is small and confounded. | Not enough for diagnostic or disease-mechanism claims. | Keep exploratory unless larger matched cohorts are added. |
| Mechanism | Feature interpretation and program enrichment exist. | Need regulatory drivers, perturbation support, and spatial localization. | Upgrade with regulon, ligand-receptor, and experimental/spatial validation. |
| Reproducibility/product | Command manifest, output provenance, artifact manifest, and manuscript rerun profile exist. | Need stable archive URIs and a clean release package. | Keep heavy artifacts out of Git; archive externally. |

## Research Milestone Dashboard

| Milestone | Status | Acceptance criterion | Required artifacts | Next action | Claim unlocked |
| --- | --- | --- | --- | --- | --- |
| R1: Regeneration-Axis Reframe | complete | ORA is framed as a healthy olfactory regenerative state axis, not an absolute age clock. | `docs/claim_ledger.md`, `docs/journal_acceptance_tracker.md`, `results/tables/ora_model_card.tsv`, `results/tables/ora_feature_interpretation.tsv`. | Maintain language discipline through final manuscript audit. | Main claim: modest, reproducible, interpretable tissue-state axis. |
| R2: Cross-Cohort Validation Atlas | in progress | Every credible public olfactory/nasal/airway aging or disease dataset is classified as direct, mapped, marker-only, context-only, or blocked; eligible datasets are processed through a shared adapter. | `results/tables/external_candidate_matrix.tsv`, `docs/journal_acceptance_tracker.md`, `docs/journal_acceptance_tracker.md`, `configs/external_datasets.yaml`, `results/tables/external_validation_evidence.tsv`. | Send/log GSE184117 label request, then promote any eligible candidate to adapter work only if the matrix class supports it. | Stronger external support if concordance survives harmonization. |
| R3: Foundation-Model Benchmark | needs work | Geneformer, scGPT, scFoundation, or CELLxGENE-derived embeddings are benchmarked against ORA/scVI/features for age-state signal and cell-state transfer. | New benchmark script, frozen model metadata, `results/tables/foundation_model_benchmark.tsv`. | Add a benchmark design doc and start with metadata-only feasibility checks. | SOTA method comparison rather than single-pipeline claim. |
| R4: Cross-Tissue Specificity | needs work | ORA features/programs are tested against non-olfactory comparators to classify signals as olfactory-specific, epithelial-regenerative, pan-airway, inflammatory, immune, or technical. | CELLxGENE Census query manifest, comparator cohort table, `results/tables/ora_cross_tissue_specificity.tsv`. | Query airway, nasal, lung, skin, gut, and immune aging resources; run feature/program overlap tests. | Specificity claim for olfactory regenerative aging if supported. |
| R5: Fate/Regeneration Dynamics | needs work | CellRank, Palantir, RNA velocity, diffusion pseudotime, or a documented no-go decision tests whether the axis aligns with basal-to-neuronal regeneration. | Lineage subset H5AD, raw-layer audit, `results/tables/regeneration_dynamics_summary.tsv`, method limitations note. | Audit raw counts/spliced layers and choose CellRank/Palantir/velocity feasibility. | Guarded regeneration trajectory support only if assumptions pass. |
| R6: Spatial/Histology Validation | blocked | Spatial transcriptomics, RNAscope, Xenium/CosMx/MERFISH, Visium, or multiplex IF localizes age-associated basal/progenitor/iOSN/sustentacular/immune programs in tissue. | Public spatial inventory or experimental design, marker panel, image/quantification table. | Search public spatial olfactory/nasal datasets; if absent, draft a validation experiment. | Tissue-localized mechanism instead of cell-suspension-only inference. |
| R7: Regulatory Mechanism Layer | needs work | Regulon/TF/program drivers are ranked and linked to age neighborhoods, lineage states, and ORA features. | SCENIC+/pySCENIC/chromVAR feasibility note, `results/tables/regulatory_driver_map.tsv`. | Run pySCENIC or a lighter TF-target enrichment prototype on lineage subsets. | Mechanistic driver hypotheses. |
| R8: Perturbation or Organoid Validation | stretch | Injury, senescence, inflammation, Notch/Wnt/EGFR, organoid, or ALI perturbation moves the ORA/regenerative axis in the predicted direction. | Perturbation inventory, experimental design, or public-data adapter. | Search public olfactory organoid/ALI/injury datasets; define minimum experiment. | Causal support if perturbation data exist or are generated. |
| R9: Clinical Phenotype Coupling | blocked | ORA or regenerative programs associate with odor function, UPSIT, smoking, COVID, CRS, biopsy site, disease duration, or other clinical covariates. | Metadata request log, covariate availability table, `results/tables/clinical_coupling_summary.tsv`. | Determine whether Gateway or external authors can provide functional/clinical covariates. | Clinical relevance beyond chronological age. |
| R10: Atlas Product Freeze | needs work | Release-ready object includes cleaned metadata, ORA score, model features, scVI/foundation embeddings, program scores, validation mappings, and a metadata dictionary. | H5AD release manifest, metadata dictionary, checksums, archive DOI/URI. | Define release schema and identify what can be redistributed. | Resource-article product claim. |

## Immediate Work Queue

| Task | Status | Artifact | Completion rule |
| --- | --- | --- | --- |
| Create SOTA research tracker | complete | `docs/sota_research_tracker.md` | Tracker exists and owns scientific upgrades. |
| Audit tracked data/database blobs | complete | Git file audit | Only `.gitkeep` placeholders are tracked under `data/` and `results/`; heavy data stay local/remote. |
| Prune stale non-pipeline docs | complete | Deleted superseded planning docs | Remaining docs are evidence, provenance, manuscript/submission assets, or generated pipeline outputs. |
| Build external candidate matrix | complete | `results/tables/external_candidate_matrix.tsv` | Each dataset has accession, assay, tissue, disease context, cohort notes, adapter status, validation class, and next action. |
| Design foundation benchmark | needs work | `docs/foundation_model_benchmark_plan.md` | Models, inputs, metrics, compute path, and failure criteria are explicit. |
| Build regeneration feature map | needs work | `results/tables/regeneration_axis_feature_map.tsv` | ORA features are assigned to basal, progenitor, iOSN, mOSN, sustentacular, immune, stress, ECM, or technical categories. |
| Audit fate-method feasibility | needs work | `results/tables/regeneration_dynamics_feasibility.tsv` | Raw counts, spliced/unspliced layers, batch structure, and lineage subsets determine allowed methods. |
| Draft spatial/perturbation validation design | needs work | `docs/spatial_perturbation_validation_plan.md` | Minimum marker panel, expected direction, samples, assay, and analysis readout are defined. |

## Claim Promotion Rules

| Claim | Main text allowed when | Otherwise |
| --- | --- | --- |
| ORA is a regenerative aging state axis | R1 remains complete and manuscript language stays modest. | Move to limitations if framed as a clock. |
| ORA is externally validated | R2 produces direct or mapped donor-level concordance with adequate metadata. | Call it small-n/contextual support only. |
| ORA is olfactory-specific | R4 shows specificity versus airway, lung, skin, gut, immune, and pan-epithelial aging comparators. | Describe shared epithelial/inflammatory aging themes. |
| Aging changes regeneration dynamics | R5 plus spatial or perturbation support shows consistent direction. | Say neighborhood or composition remodeling, not lineage flux. |
| Mechanistic TF/niche drivers are implicated | R7 links drivers to lineage states, neighborhoods, and ORA features with sensitivity checks. | Present as hypotheses. |
| Disease relevance | Larger matched disease cohorts or clinical phenotypes support it. | Keep AD/PD projection exploratory. |
| Resource-grade atlas | R10 has release schema, checksums, redistribution rules, and stable archive URI. | Keep as manuscript-internal analysis package. |

## Reviewer And PI Objection Register

| Objection | SOTA answer | Required evidence | Current status |
| --- | --- | --- | --- |
| The age model is weak. | It is not sold as an age clock; benchmark it against shuffled age, scAgeClock-style approaches, scVI, and foundation embeddings. | Model card, permutation null, calibration, foundation benchmark. | Core response ready; benchmark missing. |
| External validation is underpowered. | Exhaust public data, request labels, harmonize what is usable, and classify validation strength transparently. | External candidate matrix, adapter outputs, request log. | Needs work. |
| This is just inflammation or technical confounding. | Test cross-tissue specificity and covariate sensitivity; distinguish olfactory-regenerative, pan-epithelial, immune, and technical signatures. | Cross-tissue specificity table, matched/covariate sensitivity. | Needs work. |
| Neighborhoods are not fate. | Use fate tools only if data layers permit; otherwise add spatial/perturbation support or keep claims density-based. | Raw-layer audit, CellRank/Palantir/velocity feasibility, spatial/perturbation plan. | Needs work. |
| There is no mechanism. | Add regulon, TF-target, ligand-receptor, and perturbation-prioritized driver maps. | Regulatory driver map, niche interaction table, perturbation design. | Needs work. |
| Disease cohorts are confounded. | Do not lean on AD/PD; use disease data only as exploratory projection unless larger matched cohorts are recovered. | Disease guardrail table, label permutation, matched context. | Response ready. |
| The resource is not reusable. | Provide cleaned release object, metadata dictionary, model artifacts, checksums, and stable archive URI. | Release manifest, data dictionary, archive DOI/URI. | Needs work. |

## Data And File Hygiene Policy

| File class | Keep in Git | Keep outside Git | Clean regularly |
| --- | --- | --- | --- |
| Source code, configs, manuscript, and high-value docs | yes | no | no |
| Small final summary tables and provenance reports | yes when intentionally tracked | optional archive | no |
| Raw H5AD, processed H5AD, scVI models, large matrices, caches | no | local disk, `mia`, or durable archive | only after checksum/archive confirmation |
| Generated figures and manuscript tables | generally no unless release requires | local build outputs and archive | yes after regeneration is verified |
| Python/R/pytest/ruff caches | no | no | yes |
| Superseded planning notes | no | archive only if needed | yes |

Current tracked data audit: Git currently tracks only placeholder files under `data/` and `results/`, so the repository is not carrying large database artifacts. Heavy data cleanup should focus on ignored local files only after confirming they are archived or reproducible.

## External Sources To Track

| Source type | Why it matters | Candidate use |
| --- | --- | --- |
| Geneformer | Transformer foundation model for context-aware gene network representations. | Benchmark ORA/scVI against foundation embeddings and perturbation-relevant features. |
| scGPT | Single-cell foundation model for embedding and transfer tasks. | Compare cell-state transfer and age-state separability. |
| scFoundation | Large-scale single-cell foundation model for representation and perturbation-style tasks. | Benchmark donor/state signals and feature salience. |
| scAgeClock-style benchmarks | Single-cell aging-clock baseline. | Prevent overclaiming by comparing ORA to aging-clock expectations. |
| CellRank | Fate-mapping framework for single-cell dynamics. | Use only if input layers and assumptions pass feasibility checks. |
| Milo | Differential abundance framework for neighborhoods. | Current analysis has Milo-style and official subset sensitivity; SOTA package should preserve transparent naming. |
| CELLxGENE Census | Large public single-cell query layer. | Build cross-tissue specificity and external candidate scans. |

Reference links:

- [Geneformer](https://www.nature.com/articles/s41586-023-06139-9)
- [scGPT](https://www.nature.com/articles/s41592-024-02201-0)
- [scFoundation](https://www.nature.com/articles/s41592-024-02305-7)
- [scAgeClock](https://www.nature.com/articles/s41514-026-00379-5)
- [CellRank](https://www.nature.com/articles/s41592-021-01346-6)
- [Milo](https://www.nature.com/articles/s41587-021-01033-z)
- [CELLxGENE Census](https://chanzuckerberg.github.io/cellxgene-census/index.html)

## Next Research Decisions

| Decision | Preferred default | Why |
| --- | --- | --- |
| First SOTA analysis to implement | External candidate matrix plus regeneration feature map | Highest reviewer value and low risk. |
| First compute-heavy analysis | Foundation-model feasibility prototype | Makes the methods comparison current and credible. |
| First experimental/collaboration ask | GSE184117 labels and clinical/odor covariates | Could unlock direct validation without changing the core analysis. |
| First stretch validation design | Spatial or multiplex IF marker panel | Most direct way to prove tissue context for the regeneration axis. |
