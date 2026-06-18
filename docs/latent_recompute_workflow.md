# Latent Recompute Workflow

Updated: 2026-06-16

## Feasibility

All required latent dependencies are available.

The Gateway matrix is large, so recomputation should start with a pilot before any full-data training claim.

## Recommended Pilot Command

```bash
PYTHON=.venv/bin/python scripts/run_scvi_latent.py \
  --h5ad data/raw/gateway.h5ad \
  --out data/processed/gateway_scvi_pilot_25k.h5ad \
  --max-cells 25000 \
  --sampling-strategy stratified \
  --stratify-keys condition,fine_celltype,sex,flex_version,device_guided \
  --n-top-genes 2000 \
  --batch-key sample_id \
  --categorical-covariates flex_version,device_guided,sex \
  --hvg-flavor cell_ranger \
  --hvg-batch-key flex_version \
  --max-epochs 20 \
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
