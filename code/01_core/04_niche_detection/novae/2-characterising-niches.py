# %% Libraries
import os
import scanpy as sc
#import scvi
import numpy as np
from pathlib import Path
import pandas as pd
import anndata as ad
import rapids_singlecell as rsc
import novae
import matplotlib.pyplot as plt

# %% Setting Paths
MAIN_DIR_NAME = "cosmx_gray"
MAIN_DIR = next(p for p in Path.cwd().parents if (p / MAIN_DIR_NAME).exists()) / MAIN_DIR_NAME
os.chdir(MAIN_DIR)

# %% Setting Seed
SEED_VALUE = 42
# set NumPy RNG for consistency
np.random.seed(SEED_VALUE)

# %% object versions
CUR_OBJ_V = 'niches'
CUR_OBJ_PATH = MAIN_DIR / 'data' / 'comb' / 'h5ad' / f'comb-{CUR_OBJ_V}.h5ad'

# %% Load the unintegrated object
comb = sc.read_h5ad(CUR_OBJ_PATH)

# %% map of new names

niche_map = {
    "D1016" : "luminal",
    "D1013" : "airway_epithelium",
    "D1012" : "lymphoid_rich",
    "D1014" : "stromal-immune",
    "D1009" : "stromal-vascular",
    "D1017" : "alveolar",
}

# rename "novae_domains_7"
comb.obs['niche_names'] = comb.obs['novae_domains_6'].map(niche_map).astype('category')

# %% save annotated object
comb.write_h5ad(CUR_OBJ_PATH)