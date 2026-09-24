#!/usr/bin/env python

# %% Libraries
import os
import json
import scanpy as sc
import numpy as np
from pathlib import Path
import rapids_singlecell as rsc

# %% Main path

MAIN_DIR_NAME = "proseg_data"
MAIN_DIR = next(p for p in Path.cwd().parents if (p / MAIN_DIR_NAME).exists()) / MAIN_DIR_NAME
os.chdir(MAIN_DIR)

# %% Setting Paths -------------------------------------------------------------
CUR_OBJ_V = "integrated"
CUR_OBJ_PATH = MAIN_DIR / 'data' / 'comb' / f"h5ad/comb-{CUR_OBJ_V}.h5ad"
NEW_OBJ_V = "clustered"
NEW_OBJ_PATH = MAIN_DIR / 'data' / 'comb' / f"h5ad/comb-{NEW_OBJ_V}.h5ad"

# %% Resolution values ---------------------------------------------------------
low=0.1
high=1.5
step=0.1
RES_LIST = np.round(np.arange(low, high + 1e-9, step), 2).tolist()

# number of neighbors
N_NEIGHBORS = 15

# %% Setting Seed --------------------------------------------------------------
SEED_VALUE = 42
# set NumPy RNG for consistency
np.random.seed(42)

# %% Load comb object ----------------------------------------------------------
comb = sc.read_h5ad(CUR_OBJ_PATH)

# %% Setting constants ---------------------------------------------------------
RESOLVI_LATENT_KEY = "X_resolvi"
RESOLVI_NORMALIZED_KEY = "resolvi_normalized"
RESOLVI_NEIGHBORS_KEY = 'nb_resolvi'
RESOLVI_UMAP_KEY = 'umap_resolvi'
RESOLVI_CLUSTERS_KEY = "leiden_resolvi"

# %% Compute neighbors ---------------------------------------------------------

rsc.pp.neighbors(
    comb, 
    use_rep=RESOLVI_LATENT_KEY, 
    key_added=RESOLVI_NEIGHBORS_KEY,
    n_neighbors=N_NEIGHBORS,
    random_state=SEED_VALUE
)

# %% Compute UMAP --------------------------------------------------------------

rsc.tl.umap(
    comb, 
    neighbors_key=RESOLVI_NEIGHBORS_KEY,
    key_added=RESOLVI_UMAP_KEY,
    min_dist=0.2,
    spread=2.0,
    random_state=SEED_VALUE,
)

# Replace X_umap with integrated UMAP embeddings
comb.obsm['X_umap'] =  comb.obsm[RESOLVI_UMAP_KEY]

# %% Compute Leiden Clusters ---------------------------------------------------

rsc.tl.leiden(
    comb, 
    key_added=RESOLVI_CLUSTERS_KEY, 
    neighbors_key=RESOLVI_NEIGHBORS_KEY,
    resolution=RES_LIST, 
    n_iterations=100,
    random_state=SEED_VALUE
)

# %% Save object
comb.write_h5ad(NEW_OBJ_PATH)