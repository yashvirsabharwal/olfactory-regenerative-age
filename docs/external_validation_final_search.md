# External Validation Final Search Log

Updated: 2026-06-23

Purpose: record the final pre-submission public-data sweep for independent human olfactory/nasal validation sources. The search target was a dataset that could test donor-level ORA feature directions outside Gateway: human olfactory/nasal tissue, age or aging-relevant phenotype, expression data, donor/sample identifiers, and compatible cell labels or a realistic mapping path.

## Search Verdict

No public dataset found in this sweep is stronger than `GSE184117` for direct donor-level human olfactory aging validation. The strongest defensible position remains:

- `GSE184117` is the primary direct but small-n aging/presbyosmia validation target.
- `GSE151973` is olfactory-versus-respiratory bulk marker context only.
- `GSE290883` is a useful human olfactory-biopsy long-COVID/hyposmia single-cell context candidate, but it is not a healthy-aging replication dataset.
- `GSE290884`, `GSE235330`, `GSE324335`, `GSE151346`, and broad spatial resources remain context-only or blocked for the primary ORA claim.

## Query Families Checked

| Query family | Scope | Decision |
| --- | --- | --- |
| Human olfactory single-cell aging / presbyosmia | GEO and article/resource search for olfactory epithelial single-cell age designs. | No stronger direct dataset than `GSE184117`. |
| Human nasal / olfactory single-cell anosmia / COVID / long-COVID | Search for olfactory biopsy single-cell cohorts with disease/control design. | `GSE290883` and paired `GSE290884` added as long-COVID/hyposmia context candidates. |
| Human olfactory bulk / marker reference | Bulk olfactory-versus-respiratory references. | `GSE151973` remains marker-context only. |
| AD/PD hyposmia and neurodegeneration | Search for olfactory tissue or disease cohorts that could validate ORA disease projections. | `GSE235330` is PD/hyposmia brain-region context, not olfactory epithelium. |
| Nasal epithelial Abeta / AD-adjacent models | Search for nasal epithelial perturbation data. | `GSE324335` is immortalized nasal epithelial cell-culture context only. |
| Spatial/histology resources | HEST-style and broad spatial repositories. | No clear human olfactory/nasal spatial dataset identified for ORA feature validation. |

## Candidate Decisions

| Source | Evidence class | Inclusion decision | Reason | Next action |
| --- | --- | --- | --- | --- |
| `GSE184117` | direct, small-n; mapped/raw-adapter candidate | Keep as primary external validation target. | Human olfactory epithelium single-cell 3 prime RNA-seq with normosmic and presbyosmic adult donors. | Send/request author labels; keep scANVI and marker-reference mapping guarded. |
| `GSE151973` | marker-only/context-only | Keep as context only. | Human OE/RE bulk RNA-seq can support marker specificity but cannot test single-cell donor-level ORA features. | Use only for marker sanity language. |
| `GSE290883` | context-only; future raw-adapter candidate | Add to registry. | Human olfactory biopsies from long-COVID hyposmic and normosmic-control subjects with 5 prime GEX; disease context is not healthy aging. | Download only if the manuscript needs an additional disease/context check. |
| `GSE290884` | context-only | Add to registry. | Paired TCR-seq for the long-COVID olfactory cohort; immune repertoire context, not ORA feature validation. | Do not use for primary ORA replication. |
| `GSE235330` | context-only | Add to registry. | Human PD/hyposmia material is brain-region single-cell/snRNA-seq, not olfactory epithelium. | Use only as disease-context citation if needed. |
| `GSE324335` | context-only | Add to registry. | Immortalized human nasal epithelial cell-culture RNA-seq with SIRT3/Abeta perturbation, no donor-level olfactory-aging design. | Do not use as primary validation. |
| `GSE151346` | cross-species context-only | Keep out of human validation claims. | Mouse main olfactory epithelium single-cell context; not human donor-level validation. | Keep as optional biology context only. |
| HEST-style spatial resources | blocked/context monitor | Do not add as validation evidence yet. | Current search did not identify a clear olfactory/nasal spatial dataset with cell-state resolution. | Monitor future releases. |

## Claim Consequences

| Claim | Decision after search |
| --- | --- |
| Independent ORA replication | Not unlocked. No newly found dataset has donor-level human olfactory healthy-aging design plus labels/features. |
| GSE184117 external support | Allowed only as small-n direct/context support with marker/scANVI mapping limitations. |
| Long-COVID or anosmia biology | Context-only; can support a future discussion of olfactory injury/inflammation, not the healthy-aging axis. |
| AD/PD projection validation | Not unlocked. Disease cohorts remain exploratory and confounded. |
| Spatial localization | Defer. No ORA-ready public spatial olfactory dataset was found. |

## Sources

- `GSE184117`: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184117
- `GSE151973`: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151973
- `GSE151346`: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151346
- `GSE290883`: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE290883
- `GSE290884`: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE290884
- `GSE235330`: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235330
- `GSE324335`: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE324335
- HEST resource monitor: https://github.com/mahmoodlab/hest

