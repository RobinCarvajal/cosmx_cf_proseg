# %% Libraries
import scanpy as sc
import pandas as pd
import anndata as ad
from pathlib import Path
import os
import gc

# %% setting paths
MAIN_DIR_NAME = "proseg_data"
MAIN_DIR = next(p for p in Path.cwd().parents if (p / MAIN_DIR_NAME).exists()) / MAIN_DIR_NAME
os.chdir(MAIN_DIR)

# %% extract slide names from the slides directory
SAMPLES_DIR = Path('data/samples')

sample_names = [
    d.name for d in SAMPLES_DIR.iterdir() 
    if d.is_dir() and (
        "slide" in d.name.lower() or
        "ctrl" in d.name.lower() or
        "pwcf" in d.name.lower()
    )
]
print(f"Samples found: {sample_names}")

OBJ_VERSION = 'proseg'

# %% list to hold all sample objects
obj_list = []
for sample in sample_names:
    print(f'Loading {sample}')
    sample_path = SAMPLES_DIR / f'{sample}/h5ad/{sample}_{OBJ_VERSION}.h5ad'
    sample_obj = sc.read_h5ad(sample_path)
    # make cell_id rownames
    sample_obj.obs_names = sample_obj.obs["cell_id"].astype(str)

    obj_list.append(sample_obj)

# %% concatenate all sample objects
print("Merging sample objects...")
comb = ad.concat(
    obj_list, 
    join='outer',
)

# free memory 
del obj_list
gc.collect()

# %% save the combined object

OUTPUT_DIR = Path('data/comb/h5ad')
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
comb_path = OUTPUT_DIR / 'comb.h5ad'
comb.write_h5ad(comb_path)
# %%
