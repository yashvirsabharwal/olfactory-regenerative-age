# Spatial And Histology Validation Plan

Updated: 2026-06-25

## Decision

Direct public adult human olfactory epithelial spatial aging dataset status: `not_found`.

The best public spatial resources found are context datasets, not primary ORA validation datasets. GSE235714 is human nasal CRS/healthy-control GeoMx spatial transcriptomics; GSE292993 is human lung airway/parenchymal/vessel spatial context; GSE303809 is fetal olfactory/head-section MERFISH and is developmental only. Therefore the next strong validation step is a targeted adult olfactory histology/spatial experiment.

## Public Spatial Candidate Triage

| Dataset | Tissue | Assay | Role | Primary usable | Limitation |
| --- | --- | --- | --- | --- | --- |
| No direct adult human olfactory epithelial spatial aging dataset found | adult human olfactory epithelium | spatial transcriptomics or histology | blocking_gap | False | public data gap |
| [GSE303809 / PRJNA1297728](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE303809) | Fetal olfactory epithelium/head sections | MERFISH/spatial | developmental_spatial_context_only | False | developmental/fetal only |
| [GSE235714 / PRJNA986942](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235714) | Nasal tract / CRS tissue | NanoString GeoMx spatial | nasal_spatial_context | False | nasal/CRS context, not olfactory aging; age metadata not obvious |
| [GSE292993](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE292993) | Lung parenchyma, airway, vessel | spatial transcriptomics | airway_lung_spatial_context | False | airway/lung context, not olfactory tissue |

## Experimental Design

- Cohort: at least 10 donors per age bin if feasible: 18-40, 41-65, and 66+ years.
- Sampling: olfactory cleft/olfactory epithelium with histologic confirmation; record exact biopsy site, CRS/allergy/COVID/smoking history, sex, race/ethnicity, PMI or processing delay, and section quality.
- Exclusions or strata: active CRS/nasal polyps, acute infection, tumor-adjacent tissue, severe technical artifact, and ambiguous respiratory-only sections.
- Preferred assay: Xenium or CosMx for single-cell spatial resolution with the expanded marker panel. Practical fallback: RNAscope/IF plus a minimal marker set on adjacent sections.
- Unit of inference: donor. Regions, fields, cells, ROIs, or sections are nested observations and should not be treated as independent donors.

## Marker Panel

| Panel | Theme | Priority | Expected age direction | Markers | Readout |
| --- | --- | ---: | --- | --- | --- |
| activated_hbc_inp | activated_regeneration | 1 | context_dependent | SOX2,ASCL1,NEUROG1,NEUROD1,HES6,INSM1,MKI67,TOP2A,EGFR,KRT6A,KRT17 | Activated HBC/GBC/INP density, proliferation fraction, and distance from HBC layer. |
| hbc_quiescent_reserve | basal_reserve | 1 | positive_or_context_dependent | TP63,KRT5,KRT14,KRT15,COL17A1,ITGA6 | Basal-layer HBC area fraction and TP63/KRT5/KRT14 coexpression density. |
| immature_osn | immature_neuronal_lineage | 1 | negative | GAP43,DCX,STMN2,TUBB3,NCAM1,ELAVL4 | Immature OSN density and immature-to-mature OSN spatial ratio. |
| immune_inflammatory | immune_inflammatory | 1 | context_dependent | PTPRC,LST1,TYROBP,AIF1,C1QA,C1QB,CD3D,CD8A,NKG7,MS4A1,CXCL8,CXCL10 | Immune-cell density, epithelial infiltration, and proximity to HBC/INP and metaplastic regions. |
| mature_osn_transduction | mature_neuronal_lineage | 1 | negative | OMP,ADCY3,GNAL,GNG13,CNGA2,CNGB1,PDE1C,ANO2,RTP1,RTP2 | Mature OSN density and transduction-marker intensity in olfactory epithelial regions. |
| respiratory_metaplasia | respiratory_metaplasia | 1 | positive | FOXJ1,TPPP3,PIFO,MUC5AC,MUC5B,SPDEF,SCGB1A1,AGR2,KRT8,KRT19,CEACAM5 | Respiratory/ciliated/goblet metaplasia area fraction and boundary with olfactory regions. |
| stress_senescence_sasp | stress_senescence | 1 | positive | CDKN1A,CDKN2A,TP53,IL6,CXCL8,MMP1,MMP3,SERPINE1,CCL2,IGFBP7,DDIT3 | SASP/stress intensity by compartment and spatial overlap with metaplasia, immune, and basal regions. |
| niche_ligand_receptor | niche_signaling | 2 | context_dependent | IL17RA,IL17RC,TNF,TNFRSF1A,IFIT1,IFIT3,CXCL10,AREG,EGFR,JAG1,NOTCH1,WNT7B,FZD7,TGFB1,TGFBR2 | Ligand/receptor co-localization and sender-receiver proximity across immune, secretory, HBC, and INP regions. |
| sustentacular_barrier_detox | support_barrier_detox | 2 | context_dependent | SOX9,KRT18,CYP2A13,UGT2A1,GPX3,ALDH1A1,NQO1,HMOX1 | Sustentacular/support-cell density, detox module intensity, and epithelial-region specificity. |

## Quantitative Readouts

| Readout | Metric | Model |
| --- | --- | --- |
| compartment_area_fraction | Marker-defined compartment area fraction per section and donor. | Donor-level robust linear model against age, adjusted for sex, site, batch, and section quality when available. |
| cell_density | Marker-positive cell density per epithelial length or tissue area. | Age-bin comparison plus continuous age association; donor remains the unit of inference. |
| module_intensity | Mean marker/module intensity within annotated olfactory, respiratory, basal, neuronal, and lamina-propria regions. | Region-stratified mixed model with donor random intercept only if multiple sections per donor are present. |
| spatial_proximity | Nearest-neighbor distance or local enrichment between immune/sender cells and basal/INP/OSN receiver regions. | Permutation-calibrated spatial enrichment within sections, then donor-level age association. |
| maturation_gradient | Basal-to-apical ordering of HBC, INP, immature OSN, and mature OSN markers. | Gradient score per olfactory epithelial segment, then donor-level age association. |

## Expected Support Patterns

- Strong support: older donors show reduced immature/mature OSN spatial density or transduction intensity, altered HBC/INP localization, increased respiratory-metaplasia or stress regions, and immune/stress neighborhoods that align with ORA feature families.
- Partial support: only some ORA themes localize spatially, especially neuronal decline and respiratory metaplasia, while immune/LR signals remain context-dependent.
- Negative result: ORA-associated features do not localize to coherent epithelial or immune regions after donor-level modeling and histologic QC.

## Search Log

| Date | Resource | Query/filter | Result |
| --- | --- | --- | --- |
| 2026-06-24 | NCBI GEO DataSets / ESearch | (olfactory mucosa OR olfactory epithelium OR olfactory neuroepithelium) AND Homo sapiens AND (single cell OR single-cell OR RNA-seq OR spatial) | No direct adult human olfactory spatial aging dataset found in project registry. |
| 2026-06-24 | NCBI GEO DataSets / ESearch | (nasal OR airway OR respiratory epithelium) AND Homo sapiens AND (single cell OR single-cell OR snRNA OR spatial) AND (age OR aging OR adult) | No direct adult human olfactory spatial aging dataset found in project registry. |
| 2026-06-24 | NCBI GEO DataSets / ESearch | (olfactory OR nasal) AND Homo sapiens AND (spatial transcriptomics OR Visium OR GeoMx OR Xenium OR histology) | No direct adult human olfactory spatial aging dataset found in project registry. |
| 2026-06-25 | NCBI GEO / web search | human olfactory epithelium spatial transcriptomics; human nasal GeoMx spatial; GSE235714; GSE292993 | Confirmed nasal CRS GeoMx context dataset GSE235714 and airway/lung spatial context dataset GSE292993; no direct adult human olfactory aging spatial dataset found. |
| 2026-06-25 | Local ORA registry | public_data_exhaustion spatial/MERFISH/GeoMx candidates | GSE303809 is fetal/developmental olfactory MERFISH context only; adult olfactory spatial validation requires new data or targeted histology. |

## Claim Boundary

This plan can validate localization and compartment-level consistency of ORA biology. It cannot by itself prove lineage flux, regeneration rate, or causality; perturbation or longitudinal injury-repair data would be needed for those claims.
