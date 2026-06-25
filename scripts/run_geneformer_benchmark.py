#!/usr/bin/env python3
"""Run Geneformer donor embedding extraction for foundation-model benchmarking."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.geneformer import run_geneformer_v1_embeddings, write_geneformer_outputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", default="data/processed/foundation_benchmark_lineage_subset.h5ad")
    parser.add_argument("--model-dir", default="resources/foundation_models/geneformer_v1_10m/Geneformer-V1-10M")
    parser.add_argument("--dictionary-dir", default="resources/foundation_models/geneformer_v1_10m/geneformer/gene_dictionaries_30m")
    parser.add_argument("--max-cells", type=int, default=24_000)
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--chunk-size", type=int, default=512)
    parser.add_argument("--max-length", type=int, default=2048)
    parser.add_argument("--target-sum", type=float, default=10_000.0)
    parser.add_argument("--seed", type=int, default=20260625)
    parser.add_argument("--device", default="auto")
    parser.add_argument("--features-out", default="data/processed/geneformer_donor_features.tsv")
    parser.add_argument("--qc-out", default="results/tables/geneformer_embedding_qc.tsv")
    parser.add_argument("--coverage-out", default="results/tables/foundation_model_gene_coverage.tsv")
    parser.add_argument("--runtime-out", default="results/tables/foundation_model_runtime.tsv")
    args = parser.parse_args()

    result = run_geneformer_v1_embeddings(
        h5ad_path=args.h5ad,
        model_dir=args.model_dir,
        dictionary_dir=args.dictionary_dir,
        max_cells=args.max_cells,
        batch_size=args.batch_size,
        chunk_size=args.chunk_size,
        max_length=args.max_length,
        target_sum=args.target_sum,
        seed=args.seed,
        device=args.device,
    )
    write_geneformer_outputs(
        result,
        features_out=args.features_out,
        qc_out=args.qc_out,
        coverage_out=args.coverage_out,
        runtime_out=args.runtime_out,
    )
    qc = result.qc.iloc[0]
    print(
        "Wrote Geneformer donor features: "
        f"{args.features_out}; cells={qc['embedded_cells']}; donors={qc['donors_with_embeddings']}; "
        f"device={qc['device']}"
    )


if __name__ == "__main__":
    main()
