# %% Libraries
import os
import scanpy as sc
import numpy as np
from pathlib import Path
import rapids_singlecell as rsc

# %% Setting Paths
MAIN_DIR_NAME = "cosmx_gray"
MAIN_DIR = next(p for p in Path.cwd().parents if (p / MAIN_DIR_NAME).exists()) / MAIN_DIR_NAME
os.chdir(MAIN_DIR)

# %% Setting Seed
SEED_VALUE = 42
# set NumPy RNG for consistency
np.random.seed(SEED_VALUE)

# %% object versions
CUR_OBJ_V = 'clustered'
CUR_OBJ_PATH = MAIN_DIR / 'data' / 'comb' / 'h5ad' / f'comb-{CUR_OBJ_V}.h5ad'

# %% Load the unintegrated object
comb = sc.read_h5ad(CUR_OBJ_PATH)

# %% Map clusters to cell types (ct_resolvi)

leiden_resolvi_map = {
    # epithelial
    "8": "basal",
    "19": "secretory",
    "18": "ciliated",
    "6": "secretory-ciliated",
    "15": "AT1",
    "3": "AT2",
    # Lymphoid
    "4": "T",
    "10": "B",
    "12": "plasma",
    "14": "plasma",
    "16": "plasma",
    # Myeloid
    "0": "neutrophils",
    "2": "mast",
    "17": "alveolar_macrophages",
    "5": "macrophages",
    # Stroma
    "11": "fibroblasts",
    "7": "smooth_muscle",
    # Endothelial
    "9": "endothelial",
    # Sputum
    "13": "sputum",
    # Ambiguous
    "1": "ambiguous",
}

RESOLVI_CT_KEY = 'ct_resolvi'
RESOLVI_CLUSTERS_KEY = 'leiden_resolvi_1.0'
comb.obs['ct_resolvi'] = comb.obs[RESOLVI_CLUSTERS_KEY].map(leiden_resolvi_map).astype('category')

# %% Rename classes (LEVEL 3) - same as ct_resolvi
comb.obs['ann_lvl_3'] = comb.obs['ct_resolvi']

# %% Rename classes (LEVEL 2)

map_lvl_2 = {
    # airway epithelium 
    'basal': 'airway_epithelium',
    'secretory': 'airway_epithelium',
    'ciliated': 'airway_epithelium',
    'secretory-ciliated': 'airway_epithelium',
    # alveolar epithelium
    'AT1': 'alveolar_epithelium',
    'AT2': 'alveolar_epithelium',
    # lymphoid
    'T': 'lymphoid',
    'B': 'lymphoid',
    'plasma': 'lymphoid',
    # myeloid
    'neutrophils': 'myeloid',
    'mast': 'myeloid',
    'macrophages': 'myeloid',
    'alveolar_macrophages': 'myeloid',
    # stromal
    'fibroblasts': 'stroma',
    'smooth_muscle': 'stroma',
    # endothelial
    'endothelial': 'endothelial',
    # Sputum
    "sputum": "sputum",
    # Ambiguous
    "ambiguous": "ambiguous",
}

comb.obs['ann_lvl_2'] = comb.obs['ct_resolvi'].map(map_lvl_2).astype('category')

# %% Rename classes (LEVEL 1)

map_lvl_1 = {
    # epithelial 
    'airway_epithelium': 'epithelial',
    'alveolar_epithelium': 'epithelial',
    # immune
    'lymphoid': 'immune',
    'myeloid': 'immune',
    # stromal
    'stroma': 'stroma',
    # endothelial
    'endothelial': 'endothelial',
    # Sputum
    "sputum": "sputum",
    # Ambiguous
    "ambiguous": "ambiguous",
}

comb.obs['ann_lvl_1'] = comb.obs['ann_lvl_2'].map(map_lvl_1).astype('category')

# %% Save object
comb.write_h5ad(CUR_OBJ_PATH)
