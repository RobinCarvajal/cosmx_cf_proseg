# %% Libraries
import os
import numpy as np
import scanpy as sc
from pathlib import Path

# Note: pandas / rapids_singlecell / issparse aren't used in this block, so I've
# dropped them. Re-add any of them if a later cell needs them.

# %% Paths
# Walk up from the CWD until we find the data root (works from any nested cell folder).
PROJECT_ROOT_NAME = "cosmx_cf_proseg"
PROJECT_ROOT = next(p for p in Path.cwd().parents if (p / PROJECT_ROOT_NAME).exists()) / PROJECT_ROOT_NAME
os.chdir(PROJECT_ROOT)

# COMB -> integrated object: its .obs carries the FINAL cell-type/QC metadata
# UNINTEGRATED -> post-QC object: has the full ~6k-gene panel (the gene space we want to keep)
# OUTPUT -> object with the full gene panel and final annotations
COMB = PROJECT_ROOT / "data" / "comb" / "h5ad" / "comb-cellcharter.h5ad"
UNINTEGRATED = PROJECT_ROOT / "data" / "comb" / "h5ad" / "comb.h5ad"
OUTPUT = PROJECT_ROOT / "data" / "comb-full" / "h5ad" / "comb-full-annotated.h5ad"

# %% Reproducibility
SEED = 42
np.random.seed(SEED)

# %% Load
comb = sc.read_h5ad(COMB)     # annotations live on this object
comb_full = sc.read_h5ad(UNINTEGRATED)    # full gene panel + spatial centroids

# %% Align cells
# Keep exactly the cells present in the integrated object, in its order.
# (This also drops any unintegrated-only cells and keeps gene order intact.)
comb_full = comb_full[comb.obs_names, :].copy()

# %% Spatial coordinates
# Centroids live in comb_full's .obs and will be OVERWRITTEN when we swap in
# comb's annotation metadata below → capture them into obsm first.
comb_full.obsm["X_spatial"] = comb_full.obs[["centroid_x", "centroid_y"]].to_numpy()

# %% Attach final annotations
# comb's .obs (cell types, QC flags, …) becomes comb_full's annotation metadata.
comb_full.obs = comb.obs

# %% Normalise (in place on .X)
comb_full.layers["counts"] = comb_full.X.copy()          # raw integer counts, kept before we touch X
sc.pp.normalize_total(comb_full, target_sum=1e4)         # scale each cell to 1e4 counts
sc.pp.log1p(comb_full)                                   # log1p → .X is now the working matrix
comb_full.layers["normalized"] = comb_full.X.copy()      # copy for downstream compat; drop if unused

# %% .raw snapshot (default data source for sc.pl.*)
# Per project convention, .raw = normalised + log1p expression.
# Snapshotted AFTER the obs swap above, so comb_full.raw.obs matches the final annotations.
comb_full.raw = comb_full.copy()

# %% Write
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
comb_full.write_h5ad(OUTPUT)
print(f"Wrote {OUTPUT}  shape={comb_full.shape}  layers={list(comb_full.layers)} (+ .raw)")