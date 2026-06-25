# Perturbation, Organoid, And ALI Validation Plan

Updated: 2026-06-25

## Decision

No clean adult human olfactory aging perturbation dataset was found. GSE309325 is now adapted as olfactory-relevant organoid infection context; GSE299529 remains the strongest nasal ALI cytokine perturbation lead.

High-priority public candidates: 2. Direct adult olfactory aging perturbation datasets: 0.

## Public Candidate Triage

| Dataset | Model | Assay | Perturbation | Priority | Decision | Limitation |
| --- | --- | --- | --- | ---: | --- | --- |
| [GSE299529](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299529) | healthy donor nasal epithelial cells differentiated at air-liquid interface | scRNA-seq; processed Seurat RDS | PBS control, TNF-alpha, TGF-beta, and combined TNF-alpha/TGF-beta | 1 | high_priority_adapter_candidate | nasal respiratory model rather than olfactory epithelium; R/Seurat conversion needed locally |
| [GSE309325](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE309325) | human nasal epithelial organoids with olfactory and nasal respiratory epithelium | scRNA-seq | mock versus SARS-CoV-2 infection; 1, 2, 7, and 11 days post infection | 1 | adapter_completed_context_evidence | infection perturbation rather than aging; organoid source/developmental maturity and raw-data caveat require audit |
| [GSE286616](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE286616) | ALI-cultured human nasal epithelial cells | bulk RNA-seq | mock, rhinovirus, rhinovirus plus IRF3 inhibitor, rhinovirus plus NF-kB inhibitor, combined inhibitors | 2 | bulk_module_adapter_candidate | bulk RNA-seq; not olfactory tissue; viral-infection context |
| [GSE309353](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE309353) | well-differentiated primary nasal epithelial cells at air-liquid interface | bulk RNA-seq | mock versus RSV wild-type and NS1 mutant strains | 2 | bulk_module_adapter_candidate | viral perturbation; not olfactory or aging |
| [GSE175541](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE175541) | primary nasal epithelium co-culture with macrophages and dendritic cells | bulk RNA-seq | urban particulate matter exposure and immune co-culture across control/asthma/COPD sources | 3 | context_only_bulk_adapter_optional | asthma/COPD disease confounding; not olfactory; not aging |
| [GSE324335](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE324335) | human nasal epithelial cell culture | bulk RNA-seq | SIRT3 silencing or overexpression with amyloid-beta oligomer context | 3 | context_only_bulk_adapter_optional | immortalized/cell-culture context; not donor-level olfactory aging |
| [GSE271245](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE271245) | pediatric airway epithelial cells | bulk RNA-seq | IFN-lambda/IL-29 or IFN-beta, RSV context | 4 | low_priority_context_only | pediatric/Down syndrome airway context; not nasal/olfactory aging |

## Adapter Order

1. Audit GSE309325 first because it is the only human organoid candidate explicitly containing olfactory and respiratory epithelial compartments.
2. Build a GSE299529 adapter next if Seurat RDS conversion is practical; it is the best cytokine-driven nasal ALI single-cell perturbation candidate.
3. Use GSE286616 or GSE309353 for bulk IFN/NF-kB/viral-injury module direction only if single-cell adapters are blocked.
4. Keep GSE175541, GSE324335, and GSE271245 as lower-priority context only.

## Minimum New Experiment

| Experiment | Model | Perturbations | Timepoints | Readout | Target n |
| --- | --- | --- | --- | --- | --- |
| adult_oe_organoid_cytokine | adult donor-derived olfactory epithelial organoids or organoid-like cultures with confirmed OE and respiratory compartments | TNF-alpha, IL17A/F, IFN-beta or IFN-gamma, combined TNF/TGF-beta, and vehicle | 24h, 72h, and 7d recovery | scRNA-seq plus targeted IF/RNAscope for HBC, INP, immature OSN, mature OSN, sustentacular, respiratory, immune/stress modules | at least 6 donors, ideally balanced across young/mid/old adult age bins |
| injury_repair_axis | adult olfactory epithelial organoids or explant/ALI cultures | EGFR/AREG stimulation, Wnt activation/inhibition, Notch inhibition/activation, oxidative-stress pulse, senescence induction | baseline, injury, early recovery, late recovery | ORA/regeneration module scores, cell-state proportions, basal-to-neuronal maturation markers, and recovery trajectory | paired donor design with repeated conditions; donor remains unit of inference |

## ORA Scoring Readout

- Score ORA regeneration modules, respiratory metaplasia, IFN/TNF/IL17, senescence/SASP, oxidative stress, EGFR/AREG, Wnt, and Notch programs per condition.
- For single-cell data, summarize by donor, condition, timepoint, and harmonized epithelial state before testing perturbation effects.
- For bulk data, use module-level contrasts only; do not infer cell-state composition without single-cell or histology support.
- Treat donor or independent organoid line as the unit of inference; cells, ROIs, or technical replicates are nested observations.

## Claim Boundary

Public perturbation candidates can strengthen mechanism plausibility, especially inflammation and epithelial remodeling, but they cannot prove that ORA age associations are causal unless donor-level olfactory epithelial perturbation responses reproduce the ORA directions.

## Search Log

| Date | Resource | Query/filter | Result |
| --- | --- | --- | --- |
| 2026-06-25 | NCBI GEO DataSets / ESearch | "olfactory" AND "organoid" AND "Homo sapiens" | Returned GSE309325 and unrelated glioblastoma superseries; GSE309325 is relevant. |
| 2026-06-25 | NCBI GEO DataSets / ESearch | "olfactory epithelium" AND "injury" AND "Homo sapiens" AND "RNA-seq" | No matching GEO DataSets hit. |
| 2026-06-25 | NCBI GEO DataSets / ESearch | "nasal epithelial cells" AND "cytokine" AND "RNA-seq" | Returned GSE299529, a nasal ALI TNF/TGF-beta scRNA-seq perturbation dataset. |
| 2026-06-25 | NCBI GEO DataSets / ESearch | "nasal epithelial" AND "interferon" AND "RNA-seq" | Returned nasal/airway viral or IFN-context datasets including GSE286616, GSE271245, and GSE309353. |
| 2026-06-25 | Local ORA registry | configs/external_datasets.yaml cell-culture and nasal context entries | GSE324335 is a nasal epithelial SIRT3/amyloid-beta bulk RNA-seq context candidate. |
