#!/usr/bin/env python

'''
This script does 3 things:
1. Trains the resolvi model on the combined object
2. Extracts the latent representation and normalized expression from the model
3. Saves the updated combined object with the new data
'''


# %% Libraries
import os
import json
import scanpy as sc
import scvi
import numpy as np
from pathlib import Path
import rapids_singlecell as rsc
import pandas as pd

# %% Setting Paths
# setting paths
MAIN_DIR_NAME = "proseg_data"
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
CUR_OBJ_V = 'unintegrated'
NEW_OBJ_V = 'integrated'
CUR_OBJ_PATH = MAIN_DIR / 'data' / 'comb' / 'h5ad' / f'comb-{CUR_OBJ_V}.h5ad'
NEW_OBJ_PATH = MAIN_DIR / 'data' / 'comb' / 'h5ad' / f'comb-{NEW_OBJ_V}.h5ad'


# %% Load comb object 
comb = sc.read_h5ad(CUR_OBJ_PATH)

# %% filtering before running on resolvi
sc.pp.filter_cells(comb, min_genes=5)
sc.pp.filter_genes(comb, min_cells=3)

# %% model 

scvi.external.RESOLVI.setup_anndata(
    comb, 
    layer="counts",
    batch_key="slide_name", # huge variation at the slide level (also so that coords don't overlap)
    categorical_covariate_keys=[
        "donor",  # donor has some variation
        "run",       # run has some variation
        "modulator" # modulator might have an effect in expression
    ]
)

unsupervised_resolvi = scvi.external.RESOLVI(
    comb, 
    semisupervised=False,
    gene_likelihood="nb",
    n_layers=2,
    n_latent=10,
)

# 10 epochs should be enough for light correction
unsupervised_resolvi.train(
    max_epochs=10, # same as in default segmentation 
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


# Latent representation
latent = unsupervised_resolvi.get_latent_representation()
comb.obsm[RESOLVI_LATENT_KEY] = latent
latent.shape

# Normalized expression
norm_expr = unsupervised_resolvi.get_normalized_expression(library_size=1e4)
comb.layers[RESOLVI_NORMALIZED_KEY] = norm_expr
comb.layers[RESOLVI_NORMALIZED_KEY].shape

# %% Save integrated object as h5ad
comb.write_h5ad(NEW_OBJ_PATH)


# %% Compute neighbors on scVI latent space
'''
move this to the next script 
'''

RESOLVI_NEIGHBORS_KEY = 'nb_resolvi'
RESOLVI_UMAP_KEY = 'umap_resolvi'

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

# %% plot umap
sc.pl.umap(
    comb,
    color=['sample_name','slide_name'],
    size=0.5,
    frameon=False
)

sc.pl.umap(
    comb,
    color=['donor','run'],
    size=0.5,
    frameon=False
)
