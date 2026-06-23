# External Validation Refresh

Updated: 2026-06-23

This note records the current external-validation search state for the ORA manuscript. The main question is whether any public human olfactory/nasal dataset can test donor-level ORA feature directions independently of Gateway.

## Verdict

No newly identified dataset is currently stronger than `GSE184117` for donor-level olfactory aging validation. The final pre-submission public-data sweep is logged in `docs/external_validation_final_search.md`. The best current stance remains:

- `GSE184117`: direct but small-n presbyosmia/aging-adjacent validation target.
- `GSE151973`: bulk olfactory-versus-respiratory marker context only.
- `GSE151346`: inspected and downgraded to mouse cross-species context only; useful for olfactory epithelium cell-state/viral-entry context, not human ORA validation.
- `GSE290883`: newly logged human olfactory-biopsy long-COVID/hyposmia 5 prime GEX context candidate; not healthy-aging validation.
- `GSE290884`: paired TCR immune-context candidate for the same long-COVID olfactory cohort; not ORA feature validation.
- `GSE235330`: PD/hyposmia brain-region context only; not olfactory epithelial validation.
- `GSE324335`: nasal epithelial cell-culture Abeta/SIRT3 context only; not donor-level validation.
- Spatial transcriptomics resources such as HEST-style atlases are worth monitoring, but current search did not identify a clear human olfactory spatial dataset suitable for ORA feature validation.

## Dataset Status

| Accession / source | Assay | Tissue | Current role | ORA validation strength | Next action |
| --- | --- | --- | --- | --- | --- |
| `GSE184117` / Oliva et al. | single-cell 3' RNA-seq | olfactory epithelium | Presbyosmia and aging-adjacent validation | small-n direct/context | Keep as primary external analysis; request original cell labels if possible |
| `GSE151973` / Fodoulian et al. | bulk RNA-seq | olfactory and respiratory epithelium | OE/RE marker sanity | marker/context only | Use for olfactory-vs-respiratory marker specificity, not donor-level ORA |
| `GSE151346` / Brann et al. | mouse single-cell RNA-seq | main olfactory epithelium | COVID/anosmia mechanism context | cross-species context only | Keep out of human ORA validation claims |
| `GSE290883` / Kim et al. | single-cell 5' GEX | olfactory epithelium | Long-COVID olfactory-loss context | context/raw-adapter candidate only | Add to registry; do not promote primary ORA replication |
| `GSE290884` / Kim et al. | single-cell TCR-seq | olfactory epithelium | Long-COVID olfactory immune context | context only | Add to registry; expression validation would use `GSE290883` |
| `GSE235330` / Janssens et al. | single-cell/snRNA-seq | human brain regions | PD/hyposmia disease context | context only | Do not use as olfactory epithelial validation |
| `GSE324335` / Cartas et al. | bulk RNA-seq | immortalized nasal epithelial cells | Abeta/SIRT3 perturbation context | context only | Do not use as donor-level ORA validation |
| HEST-1k-like spatial resources | spatial transcriptomics/meta-atlas | many tissues | possible future spatial search space | not currently ORA-ready | Search metadata for olfactory/nasal tissue before investing analysis time |

## Claim Rules

- Do not call any external dataset an independent ORA replication unless it has donor/sample IDs, age or aging-relevant phenotype, expression, and Gateway-compatible cell labels or defensible reference mapping.
- COVID/anosmia and OE/RE datasets can support regeneration, sustentacular, HBC, respiratory-metaplasia, immune, and marker-specific context, but they do not validate the healthy-aging ORA axis without age design. Mouse datasets are cross-species context only.
- Spatial/histology validation remains a future opportunity. A useful public spatial dataset would need olfactory neuroepithelium or nasal olfactory-region tissue plus enough marker resolution to localize HBC, immature/mature OSN, sustentacular, secretory/glandular, and immune features.

## Sources Checked

- GSE184117 GEO record: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184117
- GSE151973 GEO record: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151973
- GSE151346 GEO record: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151346
- GSE290883 GEO record: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE290883
- GSE290884 GEO record: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE290884
- GSE235330 GEO record: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235330
- GSE324335 GEO record: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE324335
- HEST resource context: https://github.com/mahmoodlab/hest

## GSE151346 Inspection

The GEO page reports organism `Mus musculus`, and the supplementary metadata file is `GSE151346_MOE_metadata.tsv.gz`. A local inspection of that metadata found 29,585 cells plus a header row with columns including `orig_ident`, `leiden`, `leiden_name`, `UMAP_1`, and `UMAP_2`. This makes it potentially useful for cross-species olfactory epithelium cell-state context, but it is not a human validation dataset and should not be presented as evidence that ORA features replicate in human donors.
