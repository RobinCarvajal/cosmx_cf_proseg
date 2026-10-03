# %% Libraries
import os
import scanpy as sc
import scvi
import numpy as np
from pathlib import Path
import pandas as pd
import anndata as ad
import gc

# %% Setting Paths
MAIN_DIR_NAME = "proseg_data"
MAIN_DIR = next(p for p in Path.cwd().parents if (p / MAIN_DIR_NAME).exists()) / MAIN_DIR_NAME
os.chdir(MAIN_DIR)

# %% extract slide names from the slides directory
slides_dir = MAIN_DIR / 'data/slides'
slide_names = [
    d.name for d in slides_dir.iterdir() 
    if d.is_dir() and "slide" in d.name.lower()
]
print(f"Slides found: {slide_names}")

OBJ_VERSION = 'proseg'

# %% load the slide object and split it in samples

for slide in slide_names:

    # load the comb object
    slide_path = MAIN_DIR / f'data/slides/{slide}/h5ad/{slide}_{OBJ_VERSION}.h5ad'
    slide_obj = sc.read_h5ad(slide_path)

    # load the new metadata
    meta_path = MAIN_DIR / f'data/slides/{slide}/tables/{slide}_meta.tsv'
    meta = pd.read_csv(meta_path)

    # Match index type to AnnData obs_names
    meta.index = meta.index.astype(str)

    # update the metadata in the slide object
    slide_obj.obs = meta.loc[slide_obj.obs_names]
    slide_obj.obs['sample_name'] = slide_obj.obs['sample_name'].astype('category')

    # free memory
    del meta
    gc.collect()

    # extract the sample names in the slide
    sample_list = slide_obj.obs['sample_name'].cat.categories.tolist()

    for sample in sample_list:
        if sample != 'not_assigned':
            print(f'Creating object for {sample}')
            sample_obj = slide_obj[slide_obj.obs['sample_name'] == sample].copy()
            sample_path = MAIN_DIR / f'data/samples/{sample}/h5ad/{sample}_{OBJ_VERSION}.h5ad'
            sample_path.parent.mkdir(parents=True, exist_ok=True)
            sample_obj.write_h5ad(sample_path)
            print(f'Object saved at {sample_path}')

            # fee memory
            del sample_obj
            gc.collect()

    # free memory
    del slide_obj
    gc.collect()

# %%
