# External Validation Refresh

Updated: 2026-06-22

This note records the current external-validation search state for the ORA manuscript. The main question is whether any public human olfactory/nasal dataset can test donor-level ORA feature directions independently of Gateway.

## Verdict

No newly identified dataset is currently stronger than `GSE184117` for donor-level olfactory aging validation. The best current stance remains:

- `GSE184117`: direct but small-n presbyosmia/aging-adjacent validation target.
- `GSE151973`: bulk olfactory-versus-respiratory marker context only.
- `GSE151346`: newly registered single-cell olfactory mucosa/COVID-anosmia context candidate, useful for cell-state and regeneration-module sanity checks after metadata inspection.
- Spatial transcriptomics resources such as HEST-style atlases are worth monitoring, but current search did not identify a clear human olfactory spatial dataset suitable for ORA feature validation.

## Dataset Status

| Accession / source | Assay | Tissue | Current role | ORA validation strength | Next action |
| --- | --- | --- | --- | --- | --- |
| `GSE184117` / Oliva et al. | single-cell 3' RNA-seq | olfactory epithelium | Presbyosmia and aging-adjacent validation | small-n direct/context | Keep as primary external analysis; request original cell labels if possible |
| `GSE151973` / Fodoulian et al. | bulk RNA-seq | olfactory and respiratory epithelium | OE/RE marker sanity | marker/context only | Use for olfactory-vs-respiratory marker specificity, not donor-level ORA |
| `GSE151346` / Brann et al. | single-cell RNA-seq | olfactory mucosa | COVID/anosmia mechanism context | candidate marker/context | Inspect GEO supplementary files, donor/sample metadata, and annotations |
| HEST-1k-like spatial resources | spatial transcriptomics/meta-atlas | many tissues | possible future spatial search space | not currently ORA-ready | Search metadata for olfactory/nasal tissue before investing analysis time |

## Claim Rules

- Do not call any external dataset an independent ORA replication unless it has donor/sample IDs, age or aging-relevant phenotype, expression, and Gateway-compatible cell labels or defensible reference mapping.
- COVID/anosmia and OE/RE datasets can support regeneration, sustentacular, HBC, respiratory-metaplasia, immune, and marker-specific context, but they do not validate the healthy-aging ORA axis without age design.
- Spatial/histology validation remains a future opportunity. A useful public spatial dataset would need olfactory neuroepithelium or nasal olfactory-region tissue plus enough marker resolution to localize HBC, immature/mature OSN, sustentacular, secretory/glandular, and immune features.

## Sources Checked

- GSE184117 GEO record: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184117
- GSE151973 GEO record: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151973
- GSE151346 GEO record: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151346
- HEST-1k preprint/resource context: https://arxiv.org/abs/2406.16192
