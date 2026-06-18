"""Planning helpers for recomputing a non-UMAP latent space."""

from __future__ import annotations

import importlib.util
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


DEFAULT_LATENT_PACKAGES = {
    "scanpy": "scanpy",
    "scvi-tools": "scvi",
    "torch": "torch",
}
DEFAULT_SCVI_PILOT_METADATA = (
    "sample_id",
    "donor_id",
    "fine_celltype",
    "coarse_celltype",
    "flex_version",
    "device_guided",
    "sex",
    "condition",
)
DEFAULT_MARKER_PANELS = {
    "basal": ("TP63", "KRT5", "KRT14"),
    "progenitor": ("ASCL1", "NEUROG1", "NEUROD1"),
    "immature_osn": ("GAP43", "DCX", "TUBB3"),
    "mature_osn": ("OMP", "ADCY3", "GNAL"),
    "sustentacular": ("CYP2A13", "CYP2J2", "MUC1"),
    "immune": ("PTPRC", "LST1", "TYROBP"),
}


def latent_recompute_feasibility(
    schema: dict[str, Any],
    *,
    package_modules: dict[str, str] | None = None,
    n_top_genes: int = 2000,
    pilot_max_cells: int = 25_000,
) -> pd.DataFrame:
    """Return a feasibility table for scVI/scANVI latent recomputation."""

    packages = package_modules or DEFAULT_LATENT_PACKAGES
    rows = []
    missing_dependencies = False
    for package, module_name in packages.items():
        installed = importlib.util.find_spec(module_name) is not None
        missing_dependencies = missing_dependencies or not installed
        rows.append(
            {
                "check": f"dependency__{package}",
                "status": "ok" if installed else "missing",
                "detail": module_name,
                "recommendation": "available" if installed else f"Install optional latent dependency `{package}`.",
            }
        )
    n_obs = int(schema.get("n_obs") or 0)
    n_vars = int(schema.get("n_vars") or 0)
    rows.extend(
        [
            {
                "check": "input_cells",
                "status": "large",
                "detail": str(n_obs),
                "recommendation": "Use a pilot subset before full-dataset training." if n_obs > pilot_max_cells else "Pilot can use all cells.",
            },
            {
                "check": "input_genes",
                "status": "ok" if n_vars >= n_top_genes else "limited",
                "detail": str(n_vars),
                "recommendation": f"Use HVG selection with n_top_genes={min(n_top_genes, n_vars) if n_vars else n_top_genes}.",
            },
            {
                "check": "recommended_pilot",
                "status": "blocked_missing_dependencies" if missing_dependencies else "ready",
                "detail": f"max_cells={min(pilot_max_cells, n_obs) if n_obs else pilot_max_cells};n_top_genes={n_top_genes}",
                "recommendation": "Train a pilot scVI model, validate marker continuity and batch mixing, then scale up.",
            },
            {
                "check": "claim_gate",
                "status": "blocked_until_validated",
                "detail": "trajectory,milo,cnmf",
                "recommendation": "Do not promote trajectory/neighborhood/program claims until latent validation passes.",
            },
        ]
    )
    return pd.DataFrame(rows, columns=["check", "status", "detail", "recommendation"])


def render_latent_recompute_workflow(
    feasibility: pd.DataFrame,
    *,
    h5ad_path: str,
    output_h5ad: str,
    n_top_genes: int = 2000,
    pilot_max_cells: int = 250_000,
    batch_key: str = "sample_id",
    categorical_covariates: tuple[str, ...] = ("flex_version", "device_guided", "sex"),
    hvg_flavor: str = "cell_ranger",
    hvg_batch_key: str = "flex_version",
    max_epochs: int = 20,
) -> str:
    """Render a Markdown workflow for recomputing and validating scVI latent space."""

    missing = feasibility[feasibility["status"].eq("missing")]["check"].str.replace("dependency__", "", regex=False).tolist()
    dependency_line = (
        "All required latent dependencies are available."
        if not missing
        else f"Missing latent dependencies: {', '.join(missing)}."
    )
    covariates = ",".join(categorical_covariates)
    return f"""# Latent Recompute Workflow

Updated: 2026-06-16

## Feasibility

{dependency_line}

The Gateway matrix is large, so recomputation should start with a pilot before any full-data training claim.

## Recommended Pilot Command

```bash
PYTHON=.venv/bin/python scripts/run_scvi_latent.py \\
  --h5ad {h5ad_path} \\
  --out {output_h5ad} \\
  --max-cells {pilot_max_cells} \\
  --sampling-strategy stratified \\
  --stratify-keys condition,fine_celltype,sex,flex_version,device_guided \\
  --n-top-genes {n_top_genes} \\
  --batch-key {batch_key} \\
  --categorical-covariates {covariates} \\
  --hvg-flavor {hvg_flavor} \\
  --hvg-batch-key {hvg_batch_key} \\
  --max-epochs {max_epochs} \\
  --embedding-key X_scvi
```

## Validation Gates

1. Confirm `X_scvi` exists and has at least 10 dimensions.
2. Check lineage marker continuity across basal, progenitor, immature OSN, mature OSN, sustentacular, glandular, and immune states.
3. Quantify donor/sample/chemistry/collection-method mixing; reject embeddings dominated by technical axes.
4. Compare nearest-neighbor composition against known Gateway fine-cell-type labels.
5. Repeat the validation on a lineage-focused basal/progenitor/neural model.
6. Only after validation, run pseudotime, Milo, or cNMF.

## Full-Data Scaling Notes

- Use `make scvi-scaled-250k` for the first publication-scale 250k-cell, 3k-HVG stratified atlas.
- Use `make scvi-lineage-basal-neural` for the basal/progenitor/neural lineage model.
- Record runtime, peak memory, device, package versions, random seeds, and output checksums.
"""


def validate_scvi_pilot(
    h5ad_path: str | Path,
    *,
    embedding_key: str = "X_scvi",
    metadata_columns: tuple[str, ...] = DEFAULT_SCVI_PILOT_METADATA,
    min_dimensions: int = 10,
    max_validation_cells: int = 50_000,
    seed: int = 13,
) -> pd.DataFrame:
    """Validate a bounded scVI pilot H5AD without treating it as final trajectory evidence."""

    path = Path(h5ad_path)
    if not path.exists():
        return pd.DataFrame(
            [
                {
                    "check": "pilot_h5ad",
                    "status": "missing",
                    "detail": str(path),
                    "recommendation": "Run scripts/run_scvi_latent.py before validating the pilot.",
                }
            ],
            columns=["check", "status", "detail", "recommendation"],
        )
    try:
        import anndata as ad  # type: ignore
    except ModuleNotFoundError as exc:
        raise SystemExit("anndata is required to validate the scVI pilot H5AD.") from exc

    adata = ad.read_h5ad(path, backed="r")
    try:
        rows: list[dict[str, str]] = [
            {
                "check": "pilot_h5ad",
                "status": "ok",
                "detail": f"{adata.n_obs} cells x {adata.n_vars} HVGs",
                "recommendation": "Use as a bounded pilot only; full latent claims require a scaled, validated run.",
            }
        ]
        validation_indices = _validation_indices(adata.n_obs, max_cells=max_validation_cells, seed=seed)
        if embedding_key not in adata.obsm:
            rows.append(
                {
                    "check": "embedding_exists",
                    "status": "missing",
                    "detail": embedding_key,
                    "recommendation": "Rerun scVI training and confirm the embedding key is written to `.obsm`.",
                }
            )
        else:
            embedding = np.asarray(adata.obsm[embedding_key])
            n_dimensions = int(embedding.shape[1]) if embedding.ndim == 2 else 0
            finite = bool(np.isfinite(embedding).all())
            rows.extend(
                [
                    {
                        "check": "embedding_dimensions",
                        "status": "ok" if n_dimensions >= min_dimensions else "limited",
                        "detail": f"{embedding_key}:{embedding.shape}",
                        "recommendation": "Proceed to marker/mixing validation."
                        if n_dimensions >= min_dimensions
                        else f"Use at least {min_dimensions} latent dimensions for downstream validation.",
                    },
                    {
                        "check": "embedding_finite",
                        "status": "ok" if finite else "failed",
                        "detail": f"{embedding_key}:finite={finite}",
                        "recommendation": "Investigate training instability before using this embedding."
                        if not finite
                        else "No non-finite latent values detected.",
                    },
                ]
            )
            if finite and n_dimensions >= 2:
                embedding_subset = embedding[validation_indices, :]
                obs_subset = adata.obs.iloc[validation_indices].copy()
                rows.extend(_neighbor_validation_rows(embedding_subset, obs_subset))
                rows.extend(_marker_panel_rows(adata, validation_indices))
        for column in metadata_columns:
            rows.append(_metadata_validation_row(adata.obs, column))
        rows.append(
            {
                "check": "claim_gate",
                "status": "scaled_qc_seed_pending" if adata.n_obs >= 100_000 else "pilot_only",
                "detail": "trajectory,milo,cnmf",
                "recommendation": "Keep trajectory/neighborhood/program claims deferred until marker continuity and mixing diagnostics pass on a scaled latent run.",
            }
        )
    finally:
        close = getattr(adata, "file", None)
        if close is not None:
            close.close()
    return pd.DataFrame(rows, columns=["check", "status", "detail", "recommendation"])


def _validation_indices(n_obs: int, *, max_cells: int, seed: int) -> np.ndarray:
    if n_obs <= max_cells:
        return np.arange(n_obs)
    rng = np.random.default_rng(seed)
    return np.sort(rng.choice(n_obs, size=max_cells, replace=False))


def _neighbor_validation_rows(embedding: np.ndarray, obs: pd.DataFrame) -> list[dict[str, str]]:
    try:
        from sklearn.neighbors import NearestNeighbors
    except ModuleNotFoundError:
        return [
            {
                "check": "neighbor_diagnostics",
                "status": "missing_dependency",
                "detail": "sklearn",
                "recommendation": "Install scikit-learn to quantify label purity and batch mixing.",
            }
        ]
    n_neighbors = min(16, max(1, embedding.shape[0] - 2))
    if embedding.shape[0] <= 2:
        return [
            {
                "check": "neighbor_diagnostics",
                "status": "limited",
                "detail": f"cells={embedding.shape[0]}",
                "recommendation": "Use more cells for nearest-neighbor latent validation.",
            }
        ]
    neighbor_ids = NearestNeighbors(n_neighbors=n_neighbors + 1).fit(embedding).kneighbors(return_distance=False)[:, 1:]
    rows: list[dict[str, str]] = []
    for column in ("fine_celltype", "coarse_celltype"):
        if column in obs:
            purity = _mean_neighbor_same_label(obs[column].astype(str).to_numpy(), neighbor_ids)
            rows.append(
                {
                    "check": f"neighbor_label_purity__{column}",
                    "status": "ok" if purity >= 0.45 else "limited" if purity >= 0.25 else "low",
                    "detail": f"mean_same_label={purity:.3f};k={n_neighbors};cells={embedding.shape[0]}",
                    "recommendation": "Cell-state neighborhoods are coherent enough for downstream pilot interpretation."
                    if purity >= 0.45
                    else "Inspect marker continuity and label sparsity before trajectory or neighborhood claims.",
                }
            )
    for column in ("flex_version", "device_guided", "condition", "sex", "donor_id", "sample_id"):
        if column in obs:
            entropy = _mean_neighbor_entropy(obs[column].astype(str).to_numpy(), neighbor_ids)
            status = _mixing_status(column, entropy, obs[column].nunique(dropna=False))
            rows.append(
                {
                    "check": f"neighbor_mixing_entropy__{column}",
                    "status": status,
                    "detail": f"normalized_entropy={entropy:.3f};levels={obs[column].nunique(dropna=False)};k={n_neighbors}",
                    "recommendation": "No single latent-neighborhood axis should be dominated by this metadata field."
                    if status == "ok"
                    else "Review whether this metadata field is driving latent neighborhoods before making mechanistic claims.",
                }
            )
    return rows


def _marker_panel_rows(adata: Any, validation_indices: np.ndarray) -> list[dict[str, str]]:
    var_names = pd.Index(adata.var_names.astype(str))
    upper_to_name = {name.upper(): name for name in var_names}
    if "feature_name" in adata.var:
        feature_names = adata.var["feature_name"].astype(str)
        upper_to_name.update({symbol.upper(): var_name for symbol, var_name in zip(feature_names, var_names, strict=False)})
    label_column = "fine_celltype" if "fine_celltype" in adata.obs else "coarse_celltype" if "coarse_celltype" in adata.obs else ""
    rows: list[dict[str, str]] = []
    for panel, genes in DEFAULT_MARKER_PANELS.items():
        present = [upper_to_name[gene] for gene in genes if gene in upper_to_name]
        if not present:
            rows.append(
                {
                    "check": f"marker_continuity__{panel}",
                    "status": "limited",
                    "detail": "present_genes=0",
                    "recommendation": "Marker genes were not retained in the HVG subset; rerun with more HVGs or a marker-preserving feature set.",
                }
            )
            continue
        matrix = adata[validation_indices, present].X
        score = np.asarray(matrix.mean(axis=1)).ravel()
        if label_column:
            top_mask = score >= np.quantile(score, 0.9)
            labels = adata.obs.iloc[validation_indices][label_column].astype(str)
            top_labels = labels[top_mask]
            all_frequency = labels.value_counts(normalize=True)
            top_frequency = top_labels.value_counts(normalize=True)
            top_label = str(top_frequency.index[0]) if not top_frequency.empty else "none"
            enrichment = float(top_frequency.iloc[0] / max(all_frequency.get(top_label, 0.0), 1e-9)) if top_label != "none" else 0.0
            detail = f"present_genes={len(present)};top_label={top_label};top_decile_enrichment={enrichment:.2f}"
            status = "ok" if enrichment >= 2.0 else "limited"
        else:
            detail = f"present_genes={len(present)};mean_score={float(np.mean(score)):.4g}"
            status = "limited"
        rows.append(
            {
                "check": f"marker_continuity__{panel}",
                "status": status,
                "detail": detail,
                "recommendation": "Marker-high cells concentrate in an interpretable cell-state neighborhood."
                if status == "ok"
                else "Use this as a weak marker check; stronger continuity requires retained markers and label review.",
            }
        )
    return rows


def _mean_neighbor_same_label(labels: np.ndarray, neighbor_ids: np.ndarray) -> float:
    return float(np.mean(labels[neighbor_ids] == labels[:, None]))


def _mean_neighbor_entropy(labels: np.ndarray, neighbor_ids: np.ndarray) -> float:
    n_levels = len(set(labels.tolist()))
    if n_levels <= 1:
        return 0.0
    entropies = []
    denominator = np.log(n_levels)
    for row in labels[neighbor_ids]:
        _, counts = np.unique(row, return_counts=True)
        probabilities = counts / counts.sum()
        entropies.append(float(-(probabilities * np.log(probabilities)).sum() / denominator))
    return float(np.mean(entropies))


def _mixing_status(column: str, entropy: float, n_levels: int) -> str:
    if n_levels <= 1:
        return "single_level"
    if column in {"donor_id", "sample_id"} and n_levels > 100:
        return "informational" if entropy >= 0.05 else "low"
    threshold = 0.30 if column in {"flex_version", "device_guided", "condition", "sex"} else 0.20
    return "ok" if entropy >= threshold else "low"


def _metadata_validation_row(obs: pd.DataFrame, column: str) -> dict[str, str]:
    if column not in obs:
        return {
            "check": f"metadata__{column}",
            "status": "missing",
            "detail": "column_absent",
            "recommendation": "Preserve this metadata column for latent validation and manuscript reporting.",
        }
    values = obs[column].astype(str)
    counts = values.value_counts(dropna=False)
    n_levels = int(counts.shape[0])
    missing = int(values.isin(["", "nan", "None", "NA", "unknown"]).sum())
    min_cells = int(counts.min()) if not counts.empty else 0
    max_cells = int(counts.max()) if not counts.empty else 0
    if n_levels < 2:
        status = "single_level"
        recommendation = "This metadata column cannot diagnose mixing in the current pilot."
    elif min_cells < 3:
        status = "sparse_levels"
        recommendation = "Inspect rare donor/sample levels before using this pilot for neighborhood claims."
    else:
        status = "ok"
        recommendation = "Available for stratified latent validation."
    return {
        "check": f"metadata__{column}",
        "status": status,
        "detail": f"levels={n_levels};missing_or_unknown={missing};min_cells={min_cells};max_cells={max_cells}",
        "recommendation": recommendation,
    }
