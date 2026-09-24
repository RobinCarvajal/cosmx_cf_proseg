#!/usr/bin/env python

'''
This script does 3 things:
1. Trains the resolvi model on the combined object
2. Extracts the latent representation and normalized expression from the model
3. Saves the updated combined object with the new data
'''


# %% Libraries
from matplotlib import use
import os
import json
import scanpy as sc
import numpy as np
from pathlib import Path
import rapids_singlecell as rsc
import scvi
import pandas as pd

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
scvi.settings.dl_num_workers = 4
scvi.settings.dl_persistent_workers = True
scvi.settings.seed = SEED_VALUE

# %% object versions
CUR_OBJ_V = 'clustered'
NEW_OBJ_V = 'integrated'
CUR_OBJ_PATH = MAIN_DIR / 'data' / 'comb' / 'h5ad' / f'comb-{CUR_OBJ_V}.h5ad'
NEW_OBJ_PATH = MAIN_DIR / 'data' / 'comb' / 'h5ad' / f'comb-{NEW_OBJ_V}.h5ad'

# %% Results Directories
RES_DIR = "/data/cosmx_gray/code/9_extra_analysis/resolvi_de"
#RES_DIR = Path(__file__).resolve().parent
print(RES_DIR)

# %% Load comb object 
comb = sc.read_h5ad(CUR_OBJ_PATH)

# %% load resolvi model
scvi_model_path = MAIN_DIR / 'data/comb' / 'resolvi_model'
model = scvi.external.RESOLVI.load(scvi_model_path, adata=comb)

# %% model 

df = model.differential_expression(
    groupby="leiden_resolvi_1.0",
)

# save table in res directory
df.to_csv(RES_DIR / 'resolvi_de_results.csv', index=True)

# get markers 
markers = {}
for c in comb.obs['leiden_resolvi_1.0'].unique():
    cell_df = df.loc[df['group1'] == c]
    markers[c] = cell_df.index.tolist()[:10] # get top 10 markers

# dotplot 
dp = sc.pl.dotplot(
    comb, 
    markers, 
    groupby='leiden_resolvi_1.0',
    standard_scale='var',
    swap_axes=True,
    use_raw=True,
)
# save dp plot
dp.savefig(RES_DIR / 'resolvi_de_dotplot.png', dpi=300, bbox_inches='tight')