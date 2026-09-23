#!/usr/bin/env python

'''
This script does 4 things:
1. Trains the resolvi model on the combined object
2. Extracts the latent representation and normalized expression from the model
3. Computes neighbors and UMAP on the latent representation
4. Saves the updated combined object with the new data
'''


# %% Libraries
import os
import json
import scanpy as sc
import scvi
import numpy as np
from pathlib import Path
import rapids_singlecell as rsc

# %% Setting Paths
# setting paths
MAIN_DIR_NAME = "cosmx_gray"
MAIN_DIR = next(p for p in Path.cwd().parents if (p / MAIN_DIR_NAME).exists()) / MAIN_DIR_NAME
os.chdir(MAIN_DIR)

# %% Setting Seed
SEED_VALUE = 42
# set NumPy RNG for consistency
np.random.seed(42)

# %% Setting scVI parameters
scvi.settings.dl_num_workers = 64
scvi.settings.dl_persistent_workers = True
scvi.settings.seed = SEED_VALUE

# %% object versions
CUR_OBJ_V = 'v5'
#CUR_OBJ_PATH = MAIN_DIR / 'data' / 'comb' / 'h5ad' / f'comb_{CUR_OBJ_V}.h5ad'
CUR_OBJ_PATH = MAIN_DIR / 'data' / 'old_comb' / 'h5ad' / f'comb_{CUR_OBJ_V}.h5ad'
NEW_OBJ_V = 'v6'
NEW_OBJ_PATH = MAIN_DIR / 'data' / 'comb' / 'h5ad' / f'comb_{NEW_OBJ_V}.h5ad'


# %% Load comb object 
comb = sc.read_h5ad(CUR_OBJ_PATH)

# %% compute X_spatial 

comb.obsm['X_spatial'] = np.array(comb.obs[['CenterX_global_px', 'CenterY_global_px']])
# %% compute X_spatial offset

import numpy as np
import pandas as pd

coords = comb.obs[['CenterX_global_px', 'CenterY_global_px']]

sample_key = "slide_name"

# 1) Extract base coords as float arrays
x0 = coords["CenterX_global_px"].to_numpy(dtype=np.float32)
y0 = coords["CenterY_global_px"].to_numpy(dtype=np.float32)

# 2) Get sample labels as plain strings
samples = comb.obs[sample_key].astype(str)

# 3) Define a safe gap (based on global x span)
x_range = float(x0.max() - x0.min())
gap = x_range * 2   # increase to *3.0 or *5.0 if you want more separation

# 4) Stable mapping sample -> index
sample_order = pd.unique(samples)  # preserves order of appearance
sample_to_idx = {s: i for i, s in enumerate(sample_order)}

idx = samples.map(sample_to_idx).to_numpy(dtype=np.float32)

# Safety check (should never trigger)
if np.isnan(idx).any():
    bad = samples[pd.isna(samples.map(sample_to_idx))].unique()
    raise ValueError(f"Unmapped sample(s) detected: {bad}")

# 5) Create repelled coordinates
comb.obs["coords_x"] = x0 + idx * gap
comb.obs["coords_y"] = y0

# 6) Define X_spatial at the end
comb.obsm["X_spatial"] = np.array(comb.obs[["coords_x", "coords_y"]])


# save the updated object
# comb.write_h5ad(CUR_OBJ_PATH)

# %% model 

scvi.external.RESOLVI.setup_anndata(
    comb, 
    layer="counts",
    batch_key="run",
    categorical_covariate_keys=[
        "sample_name",  # sample_name has some variation
        'donor',       # donor has some variation
    ],
)

unsupervised_resolvi = scvi.external.RESOLVI(
    comb, 
    semisupervised=False,
    gene_likelihood="nb",
    n_layers=2,
    n_latent=10,
)

unsupervised_resolvi.train(
    max_epochs=10, # maybe try 100
)


# %% saving the model
models_dir = Path('data/comb')
models_dir.mkdir(parents=True, exist_ok=True)

print(models_dir)

resolvi_model_path = MAIN_DIR / 'data/comb' / 'resolvi_model'
unsupervised_resolvi.save(resolvi_model_path, overwrite=True)


# %% extract latent representation and normalized expression

# Setting constants
RESOLVI_LATENT_KEY = "X_resolvi"
RESOLVI_NORMALIZED_KEY = "resolvi_normalized"
RESOLVI_NEIGHBORS_KEY = 'nb_resolvi'
RESOLVI_UMAP_KEY = 'umap_resolvi'

# Latent representation
latent = unsupervised_resolvi.get_latent_representation()
comb.obsm[RESOLVI_LATENT_KEY] = latent
latent.shape

# Normalized expression
norm_expr = unsupervised_resolvi.get_normalized_expression(library_size=10e4)
comb.layers[RESOLVI_NORMALIZED_KEY] = norm_expr
comb.layers[RESOLVI_NORMALIZED_KEY].shape

# %% Compute neighbors on scVI latent space
rsc.pp.neighbors(
    comb,   
    use_rep=RESOLVI_LATENT_KEY, # use scVI latent space for UMAP generation
    key_added=RESOLVI_NEIGHBORS_KEY,
    n_neighbors=15,
    random_state=SEED_VALUE,
) 
# We don't use npcs. defaults to .X matrix in this case the latent the model calculated.

# %% Compute UMAP on scVI neighbors
rsc.tl.umap(
    comb, 
    neighbors_key=RESOLVI_NEIGHBORS_KEY,
    key_added=RESOLVI_UMAP_KEY,
    min_dist=0.2,
    spread=2.0,
    random_state=SEED_VALUE,
)

# %% Visualize the integrated UMAP

# Replace X_umap with integrated UMAP embeddings
comb.obsm['X_umap'] =  comb.obsm[RESOLVI_UMAP_KEY]

# %% Save integrated object as h5ad
comb.write_h5ad(NEW_OBJ_PATH)