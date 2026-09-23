# The purpose of this script is to create a minimal anndata object for cellchat 
# analysis, which will be used in the subsequent steps of the analysis pipeline. 
# The original object may contain additional data that is not necessary for 
# cellchat, so we will create a new object that only includes the essential 
# components: 
# the expression matrix (X), 
# the observation metadata (obs), 
# and the variable metadata (var). 
# This will help to reduce memory usage and improve computational efficiency 
# when running cellchat.


# %%
# %% Libraries
import os
import scanpy as sc
import numpy as np
from pathlib import Path
import anndata as ad

# %%
# %% Setting Paths
MAIN_DIR_NAME = "proseg_data"
MAIN_DIR = next(p for p in Path.cwd().parents if (p / MAIN_DIR_NAME).exists()) / MAIN_DIR_NAME
os.chdir(MAIN_DIR)

# %% Setting Seed
SEED_VALUE = 42
# set NumPy RNG for consistency
np.random.seed(SEED_VALUE)

# %% object versions
CUR_OBJ_V = 'annotated'
CUR_OBJ_PATH = MAIN_DIR / 'data' / 'comb-full' / 'h5ad' / f'comb-full-{CUR_OBJ_V}.h5ad'
CC_OBJ_PATH = MAIN_DIR / 'data' / 'comb-full' / 'h5ad' / f'comb-full-{CUR_OBJ_V}-cellchat.h5ad'

# %% Load the unintegrated object
comb = sc.read_h5ad(CUR_OBJ_PATH)

# %%
# create anndata minimal object for cellchat
cellchat_adata = ad.AnnData(
    X=comb.layers['normalized'].copy(),
    obs=comb.obs.copy(),
    var=comb.var.copy(),
)

# %%
# save 
cellchat_adata.write_h5ad(CC_OBJ_PATH)


