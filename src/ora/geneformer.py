"""Geneformer foundation-model embedding helpers."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import pickle
import time
from typing import Any

import numpy as np
import pandas as pd

from .utils import ensure_parent


@dataclass
class GeneformerEmbeddingResult:
    """Outputs from a Geneformer donor-embedding run."""

    features: pd.DataFrame
    qc: pd.DataFrame
    coverage: pd.DataFrame
    runtime: pd.DataFrame


def run_geneformer_v1_embeddings(
    *,
    h5ad_path: str | Path,
    model_dir: str | Path,
    dictionary_dir: str | Path,
    max_cells: int = 24_000,
    batch_size: int = 64,
    chunk_size: int = 512,
    max_length: int = 2048,
    target_sum: float = 10_000.0,
    seed: int = 20260625,
    device: str = "auto",
) -> GeneformerEmbeddingResult:
    """Tokenize cells with Geneformer V1 rank encoding and aggregate donor embeddings."""

    try:
        import anndata as ad  # type: ignore
        import psutil  # type: ignore
        import scipy.sparse as sp  # type: ignore
        import torch  # type: ignore
        from transformers import BertModel  # type: ignore
    except ModuleNotFoundError as exc:
        raise RuntimeError("Geneformer embedding extraction requires anndata, scipy, torch, psutil, and transformers.") from exc

    start = time.perf_counter()
    process = psutil.Process()
    token_dict, median_dict = load_geneformer_v1_dictionaries(dictionary_dir)
    adata = ad.read_h5ad(h5ad_path, backed="r")
    peak_rss = process.memory_info().rss
    try:
        gene_ids = _gene_ids(adata.var)
        gene_indices, gene_tokens, gene_medians = geneformer_gene_arrays(gene_ids, token_dict, median_dict)
        coverage = geneformer_gene_coverage(
            gene_ids,
            token_dict,
            median_dict,
            model_family="geneformer",
            checkpoint="Geneformer-V1-10M",
        )
        positions = select_geneformer_cells(adata.obs, max_cells=max_cells, seed=seed)
        obs = adata.obs.iloc[positions].copy()
        donors = sorted(obs["donor_id"].astype(str).dropna().unique().tolist())
        donor_index = {donor: idx for idx, donor in enumerate(donors)}
        dim = _hidden_size(model_dir)
        sums = np.zeros((len(donors), dim), dtype=np.float64)
        sumsq = np.zeros((len(donors), dim), dtype=np.float64)
        counts = np.zeros(len(donors), dtype=np.int64)
        model = BertModel.from_pretrained(str(model_dir), output_hidden_states=True, add_pooling_layer=False)
        run_device = resolve_device(device, torch)
        model.to(run_device).eval()

        embedded_cells = 0
        token_lengths: list[int] = []
        for start_idx in range(0, len(positions), max(1, chunk_size)):
            stop_idx = min(start_idx + max(1, chunk_size), len(positions))
            chunk_positions = positions[start_idx:stop_idx]
            chunk_obs = obs.iloc[start_idx:stop_idx]
            matrix = adata[chunk_positions, gene_indices].X
            if not sp.issparse(matrix):
                matrix = sp.csr_matrix(matrix)
            n_counts = _n_counts(chunk_obs, matrix)
            sequences, valid_local = tokenize_geneformer_matrix(
                matrix,
                n_counts=n_counts,
                gene_tokens=gene_tokens,
                gene_medians=gene_medians,
                target_sum=target_sum,
                max_length=max_length,
            )
            if not sequences:
                continue
            chunk_embeddings = embed_geneformer_sequences(
                sequences,
                model=model,
                torch_module=torch,
                device=run_device,
                batch_size=batch_size,
                pad_token_id=int(token_dict.get("<pad>", 0)),
            )
            valid_donors = chunk_obs.iloc[valid_local]["donor_id"].astype(str).to_numpy()
            donor_codes = np.array([donor_index[donor] for donor in valid_donors], dtype=np.int64)
            np.add.at(sums, donor_codes, chunk_embeddings)
            np.add.at(sumsq, donor_codes, chunk_embeddings * chunk_embeddings)
            np.add.at(counts, donor_codes, 1)
            embedded_cells += chunk_embeddings.shape[0]
            token_lengths.extend([len(seq) for seq in sequences])
            peak_rss = max(peak_rss, process.memory_info().rss)

        features = _donor_feature_table(donors, sums, sumsq, counts)
        elapsed = time.perf_counter() - start
        qc = _qc_table(
            h5ad_path=h5ad_path,
            model_dir=model_dir,
            sampled_cells=len(positions),
            embedded_cells=embedded_cells,
            donors=donors,
            counts=counts,
            token_lengths=token_lengths,
            max_length=max_length,
            batch_size=batch_size,
            device=run_device,
        )
        runtime = _runtime_table(
            status="ok",
            elapsed_seconds=elapsed,
            peak_rss_bytes=peak_rss,
            checkpoint="Geneformer-V1-10M",
            sampled_cells=len(positions),
            embedded_cells=embedded_cells,
            device=run_device,
        )
        runtime = pd.concat([runtime, deferred_foundation_model_rows()], ignore_index=True)
        return GeneformerEmbeddingResult(features=features, qc=qc, coverage=coverage, runtime=runtime)
    finally:
        close = getattr(adata, "file", None)
        if close is not None:
            close.close()


def load_geneformer_v1_dictionaries(dictionary_dir: str | Path) -> tuple[dict[str, Any], dict[str, Any]]:
    """Load Geneformer V1 token and median dictionaries."""

    root = Path(dictionary_dir)
    with (root / "token_dictionary_gc30M.pkl").open("rb") as handle:
        token_dict = pickle.load(handle)
    with (root / "gene_median_dictionary_gc30M.pkl").open("rb") as handle:
        median_dict = pickle.load(handle)
    return token_dict, median_dict


def geneformer_gene_arrays(
    gene_ids: np.ndarray,
    token_dict: dict[str, Any],
    median_dict: dict[str, Any],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return H5AD gene indices, token IDs, and median scaling factors."""

    keep = [idx for idx, gene in enumerate(gene_ids) if gene in token_dict and gene in median_dict]
    if not keep:
        raise ValueError("No H5AD genes overlap the Geneformer token and median dictionaries.")
    gene_indices = np.asarray(keep, dtype=np.int64)
    gene_tokens = np.asarray([int(token_dict[str(gene_ids[idx])]) for idx in gene_indices], dtype=np.int64)
    gene_medians = np.asarray([float(median_dict[str(gene_ids[idx])]) for idx in gene_indices], dtype=np.float32)
    gene_medians = np.where(gene_medians > 0, gene_medians, 1.0).astype(np.float32)
    return gene_indices, gene_tokens, gene_medians


def geneformer_gene_coverage(
    gene_ids: np.ndarray,
    token_dict: dict[str, Any],
    median_dict: dict[str, Any],
    *,
    model_family: str,
    checkpoint: str,
) -> pd.DataFrame:
    """Summarize Gateway gene coverage by Geneformer vocabulary."""

    genes = pd.Series(gene_ids.astype(str))
    in_token = genes.isin(token_dict.keys())
    in_median = genes.isin(median_dict.keys())
    usable = in_token & in_median
    return pd.DataFrame(
        [
            {
                "model_family": model_family,
                "checkpoint": checkpoint,
                "input_genes": int(genes.shape[0]),
                "token_vocabulary_genes": int(sum(not str(key).startswith("<") for key in token_dict)),
                "median_dictionary_genes": int(len(median_dict)),
                "genes_in_token_dictionary": int(in_token.sum()),
                "genes_in_median_dictionary": int(in_median.sum()),
                "usable_genes": int(usable.sum()),
                "usable_gene_fraction": float(usable.mean()),
                "status": "ok" if float(usable.mean()) >= 0.5 else "low_gene_coverage",
            }
        ]
    )


def select_geneformer_cells(obs: pd.DataFrame, *, max_cells: int, seed: int) -> np.ndarray:
    """Deterministically sample cells stratified by donor and fine cell type."""

    if max_cells <= 0 or max_cells >= obs.shape[0]:
        return np.arange(obs.shape[0], dtype=np.int64)
    required = {"donor_id", "fine_celltype"}
    missing = sorted(required.difference(obs.columns))
    if missing:
        raise KeyError(f"Geneformer sampling requires obs columns: {missing}")
    rng = np.random.default_rng(seed)
    work = pd.DataFrame(
        {
            "position": np.arange(obs.shape[0], dtype=np.int64),
            "donor_id": obs["donor_id"].astype(str).to_numpy(),
            "fine_celltype": obs["fine_celltype"].astype(str).to_numpy(),
        }
    )
    work["_random"] = rng.random(work.shape[0])
    n_groups = max(1, work.groupby(["donor_id", "fine_celltype"], observed=True).ngroups)
    quota = max(1, max_cells // n_groups)
    sampled = (
        work.sort_values("_random")
        .groupby(["donor_id", "fine_celltype"], observed=True, group_keys=False)
        .head(quota)
    )
    if sampled.shape[0] < max_cells:
        remaining = work.loc[~work["position"].isin(sampled["position"])]
        sampled = pd.concat([sampled, remaining.sort_values("_random").head(max_cells - sampled.shape[0])])
    elif sampled.shape[0] > max_cells:
        sampled = sampled.sort_values("_random").head(max_cells)
    return np.sort(sampled["position"].to_numpy(dtype=np.int64))


def tokenize_geneformer_matrix(
    matrix: Any,
    *,
    n_counts: np.ndarray,
    gene_tokens: np.ndarray,
    gene_medians: np.ndarray,
    target_sum: float,
    max_length: int,
) -> tuple[list[np.ndarray], list[int]]:
    """Convert a cell x gene count matrix into Geneformer V1 token sequences."""

    sequences: list[np.ndarray] = []
    valid_rows: list[int] = []
    for row_idx in range(matrix.shape[0]):
        row = matrix[row_idx].tocsr()
        if row.nnz == 0:
            continue
        denom = max(float(n_counts[row_idx]), 1.0)
        scaled = row.data.astype(np.float32) / denom * float(target_sum)
        scaled = scaled / gene_medians[row.indices]
        order = np.argsort(-scaled)[:max_length]
        tokens = gene_tokens[row.indices][order].astype(np.int64)
        if tokens.size == 0:
            continue
        sequences.append(tokens)
        valid_rows.append(row_idx)
    return sequences, valid_rows


def embed_geneformer_sequences(
    sequences: list[np.ndarray],
    *,
    model: Any,
    torch_module: Any,
    device: str,
    batch_size: int,
    pad_token_id: int,
) -> np.ndarray:
    """Run encoder-only Geneformer inference and mean-pool non-padding tokens."""

    embeddings: list[np.ndarray] = []
    for start in range(0, len(sequences), max(1, batch_size)):
        batch = sequences[start : start + max(1, batch_size)]
        max_len = max(len(seq) for seq in batch)
        input_ids = torch_module.full((len(batch), max_len), pad_token_id, dtype=torch_module.long)
        attention = torch_module.zeros((len(batch), max_len), dtype=torch_module.long)
        for idx, seq in enumerate(batch):
            input_ids[idx, : len(seq)] = torch_module.as_tensor(seq, dtype=torch_module.long)
            attention[idx, : len(seq)] = 1
        with torch_module.no_grad():
            output = model(
                input_ids=input_ids.to(device),
                attention_mask=attention.to(device),
                output_hidden_states=True,
            )
            hidden = output.hidden_states[-2]
            mask = attention.to(device).unsqueeze(-1)
            pooled = (hidden * mask).sum(dim=1) / mask.sum(dim=1).clamp_min(1)
        embeddings.append(pooled.detach().cpu().numpy().astype(np.float32))
    return np.vstack(embeddings)


def write_geneformer_outputs(
    result: GeneformerEmbeddingResult,
    *,
    features_out: str | Path,
    qc_out: str | Path,
    coverage_out: str | Path,
    runtime_out: str | Path,
) -> None:
    """Write Geneformer embedding output tables."""

    result.features.to_csv(ensure_parent(features_out), sep="\t", index=False)
    result.qc.to_csv(ensure_parent(qc_out), sep="\t", index=False)
    result.coverage.to_csv(ensure_parent(coverage_out), sep="\t", index=False)
    result.runtime.to_csv(ensure_parent(runtime_out), sep="\t", index=False)


def deferred_foundation_model_rows() -> pd.DataFrame:
    """Return explicit deferred/no-go rows for foundation models not run locally."""

    return pd.DataFrame(
        [
            {
                "model_family": "scgpt",
                "checkpoint": "whole-human",
                "status": "deferred_checkpoint_not_local",
                "elapsed_seconds": np.nan,
                "peak_rss_bytes": np.nan,
                "sampled_cells": 0,
                "embedded_cells": 0,
                "device": "",
                "notes": "Official scGPT package is PyPI-installable, but pretrained checkpoint folders are distributed through Google Drive model-zoo links; local embedding run deferred until checkpoint folder is staged.",
            },
            {
                "model_family": "scfoundation",
                "checkpoint": "100M xTrimoGene",
                "status": "deferred_weight_access_and_cli_path",
                "elapsed_seconds": np.nan,
                "peak_rss_bytes": np.nan,
                "sampled_cells": 0,
                "embedded_cells": 0,
                "device": "",
                "notes": "Official scFoundation README points to model weights/code and a newer platform/CLI path with model-license terms; local embedding run deferred until weights and CLI path are staged.",
            },
        ]
    )


def _gene_ids(var: pd.DataFrame) -> np.ndarray:
    if "ensembl_id" in var:
        return var["ensembl_id"].astype(str).to_numpy()
    return var.index.astype(str).to_numpy()


def _hidden_size(model_dir: str | Path) -> int:
    import json

    with (Path(model_dir) / "config.json").open(encoding="utf-8") as handle:
        return int(json.load(handle)["hidden_size"])


def _n_counts(obs: pd.DataFrame, matrix: Any) -> np.ndarray:
    for col in ["nCount_RNA", "n_counts", "total_counts"]:
        if col in obs:
            values = pd.to_numeric(obs[col], errors="coerce").to_numpy(dtype=np.float32)
            if np.isfinite(values).all() and np.all(values > 0):
                return values
    return np.asarray(matrix.sum(axis=1)).reshape(-1).astype(np.float32)


def _donor_feature_table(
    donors: list[str],
    sums: np.ndarray,
    sumsq: np.ndarray,
    counts: np.ndarray,
) -> pd.DataFrame:
    valid = counts > 0
    means = np.full_like(sums, np.nan, dtype=np.float64)
    sds = np.full_like(sums, np.nan, dtype=np.float64)
    means[valid] = sums[valid] / counts[valid, None]
    variance = np.maximum(sumsq[valid] / counts[valid, None] - means[valid] ** 2, 0.0)
    sds[valid] = np.sqrt(variance)
    data: dict[str, Any] = {"donor_id": donors}
    for idx in range(sums.shape[1]):
        data[f"geneformer_v1_cell_mean__dim{idx + 1:03d}"] = means[:, idx]
    for idx in range(sums.shape[1]):
        data[f"geneformer_v1_cell_sd__dim{idx + 1:03d}"] = sds[:, idx]
    return pd.DataFrame(data)


def _qc_table(
    *,
    h5ad_path: str | Path,
    model_dir: str | Path,
    sampled_cells: int,
    embedded_cells: int,
    donors: list[str],
    counts: np.ndarray,
    token_lengths: list[int],
    max_length: int,
    batch_size: int,
    device: str,
) -> pd.DataFrame:
    donor_counts = counts[counts > 0]
    lengths = np.asarray(token_lengths, dtype=float) if token_lengths else np.array([], dtype=float)
    return pd.DataFrame(
        [
            {
                "model_family": "geneformer",
                "checkpoint": "Geneformer-V1-10M",
                "h5ad_path": str(h5ad_path),
                "model_dir": str(model_dir),
                "sampled_cells": int(sampled_cells),
                "embedded_cells": int(embedded_cells),
                "n_donors": int(len(donors)),
                "donors_with_embeddings": int(np.sum(counts > 0)),
                "min_cells_per_donor": int(donor_counts.min()) if donor_counts.size else 0,
                "median_cells_per_donor": float(np.median(donor_counts)) if donor_counts.size else 0.0,
                "max_cells_per_donor": int(donor_counts.max()) if donor_counts.size else 0,
                "median_token_length": float(np.median(lengths)) if lengths.size else 0.0,
                "max_token_length": int(np.max(lengths)) if lengths.size else 0,
                "model_input_size": int(max_length),
                "batch_size": int(batch_size),
                "device": device,
            }
        ]
    )


def _runtime_table(
    *,
    status: str,
    elapsed_seconds: float,
    peak_rss_bytes: int,
    checkpoint: str,
    sampled_cells: int,
    embedded_cells: int,
    device: str,
) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "model_family": "geneformer",
                "checkpoint": checkpoint,
                "status": status,
                "elapsed_seconds": float(elapsed_seconds),
                "peak_rss_bytes": int(peak_rss_bytes),
                "sampled_cells": int(sampled_cells),
                "embedded_cells": int(embedded_cells),
                "device": device,
                "notes": "Encoder-only Geneformer V1 second-to-last-layer mean-pooled cell embeddings aggregated to donors.",
            }
        ]
    )


def resolve_device(device: str, torch_module: Any) -> str:
    """Resolve auto/cpu/mps/cuda device choice."""

    if device != "auto":
        return device
    if getattr(torch_module.backends, "mps", None) is not None and torch_module.backends.mps.is_available():
        return "mps"
    if torch_module.cuda.is_available():
        return "cuda"
    return "cpu"
