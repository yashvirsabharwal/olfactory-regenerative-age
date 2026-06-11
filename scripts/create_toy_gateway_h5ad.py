#!/usr/bin/env python3
"""Create a tiny Gateway-shaped H5AD for local smoke tests."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="data/raw/toy_gateway.h5ad")
    parser.add_argument("--donors", type=int, default=12)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    try:
        import anndata as ad  # type: ignore
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "anndata is required to create toy H5AD files. Install with `pip install -e '.[dev]'`."
        ) from exc

    rng = np.random.default_rng(args.seed)
    donors = [f"toy_donor_{i:02d}" for i in range(args.donors)]
    fine_states = [
        "Quiescent HBC",
        "Activated HBC",
        "Early INP",
        "Late INP",
        "Early iOSN",
        "Late iOSN",
        "Early mature mOSN",
        "Fully mature mOSN",
        "Stressed mOSN",
    ]
    coarse_map = {
        "Quiescent HBC": "Horizontal basal cells",
        "Activated HBC": "Horizontal basal cells",
        "Early INP": "Intermediate neuronal progenitors",
        "Late INP": "Intermediate neuronal progenitors",
        "Early iOSN": "Immature olfactory sensory neurons",
        "Late iOSN": "Immature olfactory sensory neurons",
        "Early mature mOSN": "Mature olfactory sensory neurons",
        "Fully mature mOSN": "Mature olfactory sensory neurons",
        "Stressed mOSN": "Mature olfactory sensory neurons",
    }

    obs_rows = []
    for i, donor in enumerate(donors):
        age = 35 + i * (45 / max(args.donors - 1, 1))
        disease = "healthy" if i < args.donors - 2 else ("Alzheimer's disease" if i == args.donors - 2 else "Parkinson's disease")
        chemistry = "FLEX v2" if i % 3 else "FLEX v1"
        method = "device" if i % 2 else "brush"
        for cell_idx in range(18):
            state_idx = int((cell_idx + i) % len(fine_states))
            state = fine_states[state_idx]
            obs_rows.append(
                {
                    "cell_id": f"{donor}_cell_{cell_idx:02d}",
                    "donor_id": donor,
                    "sample_id": f"toy_sample_{i:02d}",
                    "age": age,
                    "sex": "female" if i % 2 else "male",
                    "race_ethnicity": "reported",
                    "disease_condition": disease,
                    "flex_chemistry_version": chemistry,
                    "device_usage": method,
                    "coarse_cell_type": coarse_map[state],
                    "fine_cell_type": state,
                    "nCount_RNA": int(rng.integers(800, 3000)),
                    "nFeature_RNA": int(rng.integers(350, 1800)),
                    "percent.mito": float(rng.uniform(1, 10)),
                    "coarse_label_confidence": float(rng.uniform(0.75, 0.99)),
                }
            )
    obs = pd.DataFrame(obs_rows).set_index("cell_id")
    var = pd.DataFrame(index=["TP63", "ASCL1", "GAP43", "OMP", "ADCY3", "SNCA"])
    x = rng.poisson(lam=1.5, size=(obs.shape[0], var.shape[0])).astype(np.float32)
    adata = ad.AnnData(X=x, obs=obs, var=var)
    adata.obsm["X_scANVI"] = rng.normal(size=(obs.shape[0], 5)).astype(np.float32)
    output = Path(args.out)
    output.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output)
    print(f"Wrote toy H5AD: {output} ({adata.n_obs} cells x {adata.n_vars} genes)")


if __name__ == "__main__":
    main()

