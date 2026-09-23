# %% libraries 
import os
import scanpy as sc
from pathlib import Path
import numpy as np
import pandas as pd
import rapids_singlecell as rsc

# %%
MAIN_DIR_NAME = "cosmx_gray"
MAIN_DIR = next(p for p in Path.cwd().parents if (p / MAIN_DIR_NAME).exists()) / MAIN_DIR_NAME
print(MAIN_DIR)
os.chdir(MAIN_DIR)

# %% Seed value
SEED_VALUE = 42
# set NumPy RNG for consistency
np.random.seed(42)

# %% object versions
CUR_OBJ_V = 'manual-niches'
NEW_OBJ_V = 'manual-niches'
CUR_OBJ_PATH = MAIN_DIR / 'data' / 'comb' / 'h5ad' / f'comb-{CUR_OBJ_V}.h5ad'
NEW_OBJ_PATH = MAIN_DIR / 'data' / 'comb' / 'h5ad' / f'comb-{NEW_OBJ_V}.h5ad'

# %% read obj
comb = sc.read_h5ad(CUR_OBJ_PATH)

# %% Load annotations and merge with comb.obs

# Load annotations
ann_path = MAIN_DIR / "data" / "sample_sheet" / "sample_sheet.tsv"
ann = pd.read_csv(ann_path, sep="\t")

# Find columns present in both tables, except the merge key
common_cols = [
    col for col in ann.columns
    if col in comb.obs.columns and col != "sample_name"
]

# Drop overlapping columns from comb.obs
comb.obs = comb.obs.drop(columns=common_cols)

# Reset index so merge works cleanly
comb.obs = comb.obs.reset_index()

# Merge all annotation columns from ann
comb.obs = comb.obs.merge(
    ann,
    on="sample_name",
    how="left",
    validate="m:1"
)

# Restore cell_id as index
comb.obs = comb.obs.set_index("cell_id")

# %% save new obj
comb.write_h5ad(NEW_OBJ_PATH)