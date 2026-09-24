
# %% libraries 
import os
import scanpy as sc
from pathlib import Path
import numpy as np
import pandas as pd
import rapids_singlecell as rsc

# %%
MAIN_DIR_NAME = "proseg_data"
MAIN_DIR = next(p for p in Path.cwd().parents if (p / MAIN_DIR_NAME).exists()) / MAIN_DIR_NAME
print(MAIN_DIR)
os.chdir(MAIN_DIR)

# %% Seed value
SEED_VALUE = 42
# set NumPy RNG for consistency
np.random.seed(42)

# %% object versions
NEW_OBJ_V = 'unintegrated'
CUR_OBJ_PATH = MAIN_DIR / 'data' / 'comb' / 'h5ad' / 'comb.h5ad'
NEW_OBJ_PATH = MAIN_DIR / 'data' / 'comb' / 'h5ad' / f'comb-{NEW_OBJ_V}.h5ad'

# %% read obj
comb = sc.read_h5ad(CUR_OBJ_PATH)

# %% Filtering cells and genes
sc.pp.filter_cells(comb, min_genes=20)
sc.pp.filter_genes(comb, min_cells=3)

# %% Load annotations and merge with comb.obs
ann_path = MAIN_DIR / 'data' / 'sample_sheet' / 'sample_sheet.tsv'
ann = pd.read_csv(ann_path, sep='\t')

# Identify columns in ann that are NOT already in comb.obs
new_cols = [
    c for c in ann.columns 
    if c not in comb.obs.columns or 
    c == 'sample_name'
]

# reset the index of comb.obs to enable merging
comb.obs.reset_index(drop=True, inplace=True)

# Merge using only non-duplicate columns
comb.obs = comb.obs.merge(
    ann[new_cols],
    on='sample_name',
    how='left',
    validate='m:1'
)

# set index to cell_id if not already set
if comb.obs.index.name != 'cell_id':
    comb.obs.set_index('cell_id', inplace=True)

# %% define X_spatial
comb.obsm["X_spatial"] = np.array(
    comb.obs[['centroid_x', 'centroid_y']]
)

# %% Normalization and log1p transformation
comb.layers['counts'] = comb.X.copy() # saving the raw counts in a layer
sc.pp.normalize_total(comb, target_sum=1e4)
sc.pp.log1p(comb)
comb.raw = comb # this is just a snapshot of normalized+1logp counts. Used by default for plots.

# %% Highly variable genes
sc.pp.highly_variable_genes(
    comb,
    n_top_genes=2000, # change to 2000 and re-train the model
    subset=True, # do not subset to highly variable genes
    layer="counts",
    flavor="seurat_v3",
    batch_key="run",
)

sc.pl.highly_variable_genes(comb)

# %% PCA calculation
rsc.pp.pca(comb)
sc.pl.pca_variance_ratio(comb, n_pcs=50, log=True)

# %% Neighbors calculation 
UNINTEGRATED_NEIGHBORS_KEY = 'nb_unint'
rsc.pp.neighbors(
    comb,
    n_pcs=30,
    key_added=UNINTEGRATED_NEIGHBORS_KEY,
    random_state=SEED_VALUE,
)

comb.uns.keys() # check the keys in uns

# %% UMAP calculation 
UNINTEGRATED_UMAP_KEY = 'umap_unint'
rsc.tl.umap(
    comb, 
    key_added=UNINTEGRATED_UMAP_KEY, 
    neighbors_key=UNINTEGRATED_NEIGHBORS_KEY,
    random_state=SEED_VALUE,
)

comb.uns.keys() # check the keys in uns

# %% create X_umap for sc.pl.umap to work and assign unintegrated embeddings
comb.obsm['X_umap'] = comb.obsm[UNINTEGRATED_UMAP_KEY]

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

# %% Check object before saving
comb

# %% Save unintegrated object as h5ad
comb.write_h5ad(NEW_OBJ_PATH)
