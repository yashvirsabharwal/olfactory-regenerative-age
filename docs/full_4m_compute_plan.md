# Full 4M scVI Compute Plan

Updated: 2026-06-18

## Decision

Use the remote server path for full 4M-cell scVI work. The direct raw-H5AD backed row/column slicing route exceeded practical memory, but the chunked reduced-H5AD route completed successfully on `mia`: all 4,028,275 cells were reduced to a 3,003-gene HVG/marker feature set and scVI was trained on all cells.

## Current SSH Status

`mia.ninds.nih.gov` is reachable. Password SSH should be handled through an SSH control socket so the password is entered once in a normal terminal and subsequent non-interactive `ssh`/`rsync` calls reuse the authenticated connection.

If the control socket is not already active, establish one of:

- a valid SSH key accepted by `sabharwaly2@mia.ninds.nih.gov`;
- a Kerberos/GSSAPI session usable by SSH from this workstation;
- or an interactive login path, followed by running the commands manually on the server.

For password-based SSH, use:

```bash
mkdir -p ~/.ssh/cm
ssh -MNf \
  -o ControlMaster=yes \
  -o ControlPath=~/.ssh/cm/%r@%h:%p \
  -o ControlPersist=8h \
  sabharwaly2@mia.ninds.nih.gov
```

Enter the password at the prompt. Then run the remote helper with the same control path:

```bash
SSH_OPTS="-o ControlPath=~/.ssh/cm/%r@%h:%p" scripts/remote_full_scvi.sh check
SSH_OPTS="-o ControlPath=~/.ssh/cm/%r@%h:%p" scripts/remote_full_scvi.sh push
SSH_OPTS="-o ControlPath=~/.ssh/cm/%r@%h:%p" scripts/remote_full_scvi.sh launch
```

Close the master connection when done:

```bash
ssh -O exit -o ControlPath=~/.ssh/cm/%r@%h:%p sabharwaly2@mia.ninds.nih.gov
```

## Full-Model Targets

The repo now has first-class targets:

```bash
PYTHON=.venv/bin/python make scvi-reduced-4m
PYTHON=.venv/bin/python make scvi-full-4m-reduced
PYTHON=.venv/bin/python make scvi-full-4m
PYTHON=.venv/bin/python make scvi-full-validation
PYTHON=.venv/bin/python make output-provenance
```

Expected outputs:

- `data/processed/gateway_hvg3003_4m.h5ad`
- `data/processed/gateway_scvi_full_4m_reduced.h5ad`
- `results/models/gateway_scvi_full_4m_reduced/`
- `results/tables/scvi_full_4m_reduced_validation.tsv`

The preferred all-cell route on `mia` is now:

1. `scvi-reduced-4m`: read `data/raw/gateway.h5ad` in small contiguous row chunks, keep the 3,003-gene HVG/marker set, write temporary reduced chunk H5ADs, and concatenate them on disk into `data/processed/gateway_hvg3003_4m.h5ad`.
2. `scvi-full-4m-reduced`: train scVI on all cells from the reduced all-cell object.

This preserves all 4,028,275 cells while avoiding the memory-hostile backed row/column fancy-index operation that blocked direct training from the raw H5AD.

If the 3,003-gene all-cell run exceeds practical memory during H5AD materialization, use the memory-constrained fallback:

```bash
PYTHON=.venv/bin/python make scvi-full-4m-safe
```

This still uses all Gateway cells, but limits the model to a marker-preserving 1,500-gene feature set ranked from the validated 250k scVI reference by marker status, HVG batch support, and normalized dispersion. Its expected outputs are:

- `data/processed/gateway_scvi_full_4m_safe.h5ad`
- `results/models/gateway_scvi_full_4m_safe/`

If the all-cell fallback still exceeds memory, run the large stratified atlas:

```bash
PYTHON=.venv/bin/python make scvi-scaled-1m
```

This uses 1,000,000 stratified cells and the full 3,003-gene selected feature set. It is the preferred defensible sketch when the current AnnData/scVI path cannot materialize all 4,028,275 cells on the available no-swap server.

If the 1M sampled materialization also approaches the memory ceiling, use:

```bash
PYTHON=.venv/bin/python make scvi-scaled-500k
```

This preserves the full selected 3,003-gene feature set while halving the sampled-cell memory pressure.

If 500k also approaches the memory ceiling, run a second 250k seed-stability model instead:

```bash
PYTHON=.venv/bin/python make scvi-scaled-250k-seed23
```

This does not replace an eventual full atlas, but it gives a publication-useful stability check at the largest scale already known to complete with the current in-memory materialization route.
This target intentionally follows the known-good row-sampling plus in-memory HVG-selection route, rather than preselecting columns from the backed CSR H5AD, because backed row/column fancy indexing showed high transient memory on `mia`.

## Remote Runner

After SSH works:

```bash
scripts/remote_full_scvi.sh check
scripts/remote_full_scvi.sh push
scripts/remote_full_scvi.sh launch
scripts/remote_full_scvi.sh status
scripts/remote_full_scvi.sh fetch
```

Optional overrides:

```bash
REMOTE=sabharwaly2@mia.ninds.nih.gov \
REMOTE_DIR=/path/to/large/scratch/olfactory \
PYTHON_BIN=python3 \
scripts/remote_full_scvi.sh push
```

Use a large scratch filesystem if home storage is quota-limited. The `push` action syncs source code and `data/raw/gateway.h5ad`; the `fetch` action retrieves only the full-model H5AD, model directory, validation table, provenance table, and log.

## Acceptance Checks

- Full-model H5AD has finite `X_scvi` with at least 10 dimensions. Completed for `gateway_scvi_full_4m_reduced.h5ad`, `X_scvi:(4028275, 10)`.
- `results/tables/scvi_full_4m_reduced_validation.tsv` exists and passes the same gates used for 250k validation.
- Output provenance reports 0 missing outputs after fetch.
- GSE184117 is remapped against the full reference only after full-model validation passes.

## Fallback

The first all-cell 3,003-gene launch on `mia` reached roughly 76 GB RSS during backed H5AD materialization, exited before GPU training, and produced no full-model output. The 1,500-gene all-cell fallback reached roughly 103 GB RSS while still materializing and was stopped before exhausting the 125 GB no-swap server. The 1M-cell and 500k-cell stratified runs both reached roughly 106 GB RSS during sampled materialization and were also stopped. The successful route was the chunked reduced-H5AD workflow: read contiguous raw rows, subset genes inside each small in-memory chunk, concatenate chunks on disk after installing `dask[array]`, then train scVI on the reduced all-cell object.
