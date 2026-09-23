
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

# %% Load comb object 

CUR_OBJ_V = 'v4'
CUR_OBJ_PATH = MAIN_DIR / 'data' / 'comb' / 'h5ad' / f'comb_{CUR_OBJ_V}.h5ad'
NEW_OBJ_V = 'v5'
NEW_OBJ_PATH = MAIN_DIR / 'data' / 'comb' / 'h5ad' / f'comb_{NEW_OBJ_V}.h5ad'

comb = sc.read_h5ad(CUR_OBJ_PATH)

# %% Load scanvi model
scanvi_model_path = MAIN_DIR / 'data/comb' / 'scanvi_model'
scanvi_model = scvi.model.SCANVI.load(
    dir_path=scanvi_model_path, 
    adata=comb,
    accelerator='gpu',
)

# %% Get latent representation
SCANVI_LATENT_KEY = "X_scANVI"
comb.obsm[SCANVI_LATENT_KEY] = scanvi_model.get_latent_representation(comb)

# %% Setting key names for neighbors and umap
SCANVI_NEIGHBORS_KEY = "nb_scANVI"
SCANVI_UMAP_KEY = "umap_scANVI"

# %% Compute neighbors
# use scANVI latent space for Neighbors calculation
rsc.pp.neighbors(
    comb, 
    use_rep=SCANVI_LATENT_KEY, 
    key_added=SCANVI_NEIGHBORS_KEY,
    n_neighbors=15,
    random_state=SEED_VALUE,
)

# %% Predicting cell types
SCANVI_PRED_KEY = "ct_scANVI"
preds = scanvi_model.predict(comb)      # predicted labels
comb.obs[SCANVI_PRED_KEY] = preds

# %% Run UMAP
rsc.tl.umap(
    comb,
    neighbors_key=SCANVI_NEIGHBORS_KEY,
    key_added=SCANVI_UMAP_KEY,
    min_dist=0.2,
    spread=2.0,
    random_state=SEED_VALUE,
)

# set the default umap to the scanvi umap
comb.obsm['X_umap'] =  comb.obsm[SCANVI_UMAP_KEY]

# %% plot UMAP
sc.pl.umap(
    comb,
    color=['ct_scVI','ct_scANVI','leiden_scVI'],
    frameon=False,
    ncols=1,
    legend_loc='on data'
)

# %% check marker expression

# load all_markers from json file
all_markers_path = 'code/canonical_markers/all_markers.json'
with open(all_markers_path, 'r') as f:
    all_markers = json.load(f)

# dotplot
dp = sc.pl.dotplot(comb, all_markers, groupby='ct_scANVI', standard_scale='var', return_fig=True)
dp.add_totals().show()

# %% save integrate object as h5ad
comb.write_h5ad(MAIN_DIR / 'data/comb/h5ad' / 'comb_v5.h5ad')

# %% save metadata
meta = comb.obs.copy()
meta.to_csv(MAIN_DIR / 'data/comb/meta' / 'comb-meta.tsv', sep='\t')

# %%
