#!/usr/bin/env python3
"""Train a Gateway scANVI reference and map GSE184117 query cells."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.reference_mapping import mapped_label_donor_features, mapping_qc_by_sample, normalized_entropy
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-h5ad", default="data/processed/gateway_scvi_stratified_250k.h5ad")
    parser.add_argument("--query-h5ad", default="data/processed/gse184117_marker_mapped.h5ad")
    parser.add_argument("--model-dir", default="results/models/gateway_scanvi_reference")
    parser.add_argument("--query-out", default="data/processed/gse184117_scanvi_mapped.h5ad")
    parser.add_argument("--qc-out", default="results/tables/external_scanvi_mapping_qc.tsv")
    parser.add_argument("--donor-features-out", default="data/processed/gse184117_scanvi_donor_features.tsv")
    parser.add_argument("--metadata-out", default="results/tables/gateway_scanvi_reference_metadata.tsv")
    parser.add_argument("--label-key", default="fine_celltype")
    parser.add_argument("--batch-key", default="sample_id")
    parser.add_argument("--categorical-covariates", default="flex_version,device_guided,sex")
    parser.add_argument("--unlabeled-category", default="Unknown")
    parser.add_argument("--max-epochs-scvi", type=int, default=30)
    parser.add_argument("--max-epochs-scanvi", type=int, default=20)
    parser.add_argument("--max-epochs-query", type=int, default=10)
    parser.add_argument("--seed", type=int, default=23)
    parser.add_argument("--overwrite-model", action="store_true")
    parser.add_argument("--reuse-reference", action="store_true")
    args = parser.parse_args()

    try:
        import anndata as ad  # type: ignore
        import scvi  # type: ignore
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "anndata and scvi-tools are required for scANVI/scArches mapping. "
            "Install `python -m pip install -e .[latent]` first."
        ) from exc

    scvi.settings.seed = args.seed
    model_dir = Path(args.model_dir)
    reference_adata = None
    if args.reuse_reference and model_dir.exists():
        reference_adata = ad.read_h5ad(args.reference_h5ad)
        scanvi_model = scvi.model.SCANVI.load(str(model_dir), adata=reference_adata)
    else:
        reference_adata = ad.read_h5ad(args.reference_h5ad)
        _prepare_reference_obs(
            reference_adata.obs,
            label_key=args.label_key,
            unlabeled_category=args.unlabeled_category,
        )
        batch_key = args.batch_key if args.batch_key in reference_adata.obs else None
        covariates = [
            key
            for key in _split_arg(args.categorical_covariates)
            if key in reference_adata.obs and key != batch_key and key != args.label_key
        ]
        scvi.model.SCVI.setup_anndata(
            reference_adata,
            batch_key=batch_key,
            labels_key=args.label_key,
            categorical_covariate_keys=covariates or None,
        )
        scvi_model = scvi.model.SCVI(reference_adata)
        scvi_model.train(max_epochs=args.max_epochs_scvi, early_stopping=True)
        scanvi_model = scvi.model.SCANVI.from_scvi_model(
            scvi_model,
            labels_key=args.label_key,
            unlabeled_category=args.unlabeled_category,
        )
        scanvi_model.train(max_epochs=args.max_epochs_scanvi, early_stopping=True)
        model_dir.mkdir(parents=True, exist_ok=True)
        scanvi_model.save(str(model_dir), overwrite=args.overwrite_model, save_anndata=False)

    query = ad.read_h5ad(args.query_h5ad)
    _prepare_query_obs(
        query.obs,
        label_key=args.label_key,
        unlabeled_category=args.unlabeled_category,
        categorical_covariates=_split_arg(args.categorical_covariates),
    )
    query_model = scvi.model.SCANVI.load_query_data(
        query,
        reference_model=scanvi_model,
        inplace_subset_query_vars=True,
    )
    query_model.train(max_epochs=args.max_epochs_query, early_stopping=True)
    query_adata = query_model.adata
    query_adata.obsm["X_scanvi"] = query_model.get_latent_representation()
    predictions = _as_prediction_array(query_model.predict())
    probabilities = query_model.predict(soft=True)
    probability_frame = _as_probability_frame(probabilities)
    query_adata.obs["scanvi_predicted_label"] = np.asarray(predictions).astype(str)
    query_adata.obs["scanvi_label_confidence"] = probability_frame.max(axis=1).to_numpy(dtype=float)
    query_adata.obs["scanvi_label_entropy"] = normalized_entropy(probability_frame)
    for label in probability_frame.columns:
        query_adata.obs[f"scanvi_prob__{_safe_label(label)}"] = probability_frame[label].to_numpy(dtype=float)

    qc = mapping_qc_by_sample(
        query_adata.obs,
        label_column="scanvi_predicted_label",
        confidence_column="scanvi_label_confidence",
        entropy_column="scanvi_label_entropy",
    )
    donor_features = mapped_label_donor_features(
        query_adata.obs,
        label_column="scanvi_predicted_label",
        confidence_column="scanvi_label_confidence",
    )
    metadata = _reference_metadata(
        reference_adata,
        query_adata,
        args=args,
        model_dir=model_dir,
        label_probabilities=probability_frame,
    )
    ensure_parent(args.query_out)
    query_adata.write_h5ad(args.query_out)
    qc.to_csv(ensure_parent(args.qc_out), sep="\t", index=False)
    donor_features.to_csv(ensure_parent(args.donor_features_out), sep="\t", index=False)
    metadata.to_csv(ensure_parent(args.metadata_out), sep="\t", index=False)
    print(f"Wrote scANVI query mapping: {args.query_out} ({query_adata.n_obs} cells x {query_adata.n_vars} genes)")
    print(f"Wrote scANVI mapping QC: {args.qc_out} ({qc.shape[0]} samples)")
    print(f"Wrote scANVI donor features: {args.donor_features_out} ({donor_features.shape[0]} rows)")
    print(f"Wrote scANVI reference metadata: {args.metadata_out}")


def _prepare_reference_obs(obs: pd.DataFrame, *, label_key: str, unlabeled_category: str) -> None:
    if label_key not in obs:
        raise ValueError(f"Reference AnnData is missing label column {label_key!r}.")
    labels = obs[label_key].astype(str).replace({"": unlabeled_category, "nan": unlabeled_category, "None": unlabeled_category})
    obs[label_key] = pd.Categorical(labels)
    if unlabeled_category not in obs[label_key].cat.categories:
        obs[label_key] = obs[label_key].cat.add_categories([unlabeled_category])


def _prepare_query_obs(
    obs: pd.DataFrame,
    *,
    label_key: str,
    unlabeled_category: str,
    categorical_covariates: tuple[str, ...],
) -> None:
    obs[label_key] = pd.Categorical([unlabeled_category] * obs.shape[0])
    for column in categorical_covariates:
        if column not in obs:
            obs[column] = "query_unknown"
        obs[column] = obs[column].astype(str).replace({"": "query_unknown", "nan": "query_unknown", "None": "query_unknown"})


def _as_probability_frame(probabilities) -> pd.DataFrame:
    if isinstance(probabilities, tuple):
        probabilities = probabilities[0]
    if isinstance(probabilities, pd.DataFrame):
        return probabilities.reset_index(drop=True)
    array = np.asarray(probabilities, dtype=float)
    return pd.DataFrame(array, columns=[f"label_{idx}" for idx in range(array.shape[1])])


def _as_prediction_array(predictions) -> np.ndarray:
    if isinstance(predictions, tuple):
        predictions = predictions[0]
    if isinstance(predictions, pd.Series):
        return predictions.astype(str).to_numpy()
    if isinstance(predictions, pd.DataFrame):
        if predictions.shape[1] != 1:
            raise ValueError(f"Expected one prediction column, found {predictions.shape[1]}.")
        return predictions.iloc[:, 0].astype(str).to_numpy()
    return np.asarray(predictions).astype(str)


def _reference_metadata(
    reference_adata,
    query_adata,
    *,
    args: argparse.Namespace,
    model_dir: Path,
    label_probabilities: pd.DataFrame,
) -> pd.DataFrame:
    labels = reference_adata.obs[args.label_key].astype(str) if reference_adata is not None and args.label_key in reference_adata.obs else pd.Series(dtype=str)
    payload = {
        "reference_h5ad": args.reference_h5ad,
        "query_h5ad": args.query_h5ad,
        "model_dir": str(model_dir),
        "reference_cells": int(reference_adata.n_obs) if reference_adata is not None else 0,
        "reference_genes": int(reference_adata.n_vars) if reference_adata is not None else 0,
        "query_cells": int(query_adata.n_obs),
        "query_genes_after_intersection": int(query_adata.n_vars),
        "label_key": args.label_key,
        "n_reference_labels": int(labels.nunique()) if not labels.empty else 0,
        "n_probability_labels": int(label_probabilities.shape[1]),
        "batch_key": args.batch_key,
        "categorical_covariates": args.categorical_covariates,
        "max_epochs_scvi": int(args.max_epochs_scvi),
        "max_epochs_scanvi": int(args.max_epochs_scanvi),
        "max_epochs_query": int(args.max_epochs_query),
        "seed": int(args.seed),
        "mapping_method": "scANVI_load_query_data_scArches_style",
    }
    payload["label_counts_json"] = json.dumps(labels.value_counts().to_dict(), sort_keys=True) if not labels.empty else "{}"
    return pd.DataFrame([payload])


def _safe_label(label: str) -> str:
    return str(label).strip().lower().replace(" ", "_").replace("-", "_").replace("/", "_")


def _split_arg(value: str) -> tuple[str, ...]:
    return tuple(part.strip() for part in value.split(",") if part.strip())


if __name__ == "__main__":
    main()
