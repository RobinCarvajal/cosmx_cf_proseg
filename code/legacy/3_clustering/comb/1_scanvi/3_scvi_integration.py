
# %% Libraries
import os
import scanpy as sc
import scvi
import numpy as np
from pathlib import Path
import pandas as pd
import anndata as ad
import rapids_singlecell as rsc

# %% Setting Paths
MAIN_DIR_NAME = "cosmx_gray"
MAIN_DIR = next(p for p in Path.cwd().parents if (p / MAIN_DIR_NAME).exists()) / MAIN_DIR_NAME
os.chdir(MAIN_DIR)

# %% Setting Seed
SEED_VALUE = 42
# set NumPy RNG for consistency
np.random.seed(42)
scvi.settings.seed = SEED_VALUE

# %% object versions
CUR_OBJ_V = 'v1'
CUR_OBJ_PATH = MAIN_DIR / 'data' / 'comb' / 'h5ad' / f'comb_{CUR_OBJ_V}.h5ad'
NEW_OBJ_V = 'v2'
NEW_OBJ_PATH = MAIN_DIR / 'data' / 'comb' / 'h5ad' / f'comb_{NEW_OBJ_V}.h5ad'

# %% Load the unintegrated object
comb = sc.read_h5ad(CUR_OBJ_PATH)

# %% Load the model
scvi_model_path = MAIN_DIR / 'data/comb' / 'scvi_model'
model = scvi.model.SCVI.load(scvi_model_path, adata=comb)

# %% Setting constants
SCVI_LATENT_KEY = "X_scVI"
SCVI_NORMALIZED_KEY = "scVI_normalized"
SCVI_NEIGHBORS_KEY = 'nb_scVI'
SCVI_UMAP_KEY = 'umap_scVI'

# %% Get the model outputs

# Latent representation
latent = model.get_latent_representation()
comb.obsm[SCVI_LATENT_KEY] = latent
latent.shape

# Normalized expression
comb.layers[SCVI_NORMALIZED_KEY] = model.get_normalized_expression(library_size=10e4)
comb.layers[SCVI_NORMALIZED_KEY].shape

# %% Compute neighbors on scVI latent space
rsc.pp.neighbors(
    comb,   
    use_rep=SCVI_LATENT_KEY, # use scVI latent space for UMAP generation
    key_added=SCVI_NEIGHBORS_KEY,
    n_neighbors=15,
    random_state=SEED_VALUE,
) 
# We don't use npcs. defaults to .X matrix in this case the latent the model calculated.

# %% Compute UMAP on scVI neighbors
rsc.tl.umap(
    comb, 
    neighbors_key=SCVI_NEIGHBORS_KEY,
    key_added=SCVI_UMAP_KEY,
    min_dist=0.2,
    spread=2.0,
    random_state=SEED_VALUE,
)

# %% Visualize the integrated UMAP

# Replace X_umap with integrated UMAP embeddings
comb.obsm['X_umap'] =  comb.obsm[SCVI_UMAP_KEY]

sc.pl.umap(
    comb,
    color=["sample_name", "slide_name"],
    size=0.5,
    frameon=False
)

sc.pl.umap(
    comb,
    color=["donor", "run"],
    size=0.5,
    frameon=False
)

# %% Save integrated object as h5ad
comb.write_h5ad(NEW_OBJ_PATH)

