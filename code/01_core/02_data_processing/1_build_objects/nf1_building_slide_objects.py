# %% load libraries
from pathlib import Path
import os
import scanpy as sc
import squidpy as sq
import gc

import argparse
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument(
    "--input", 
    required=True, 
    type=Path,
    help="Path to the directory containing the slide flat files")
p.add_argument(
    "--output", 
    required=True, 
    type=Path,
    help="Path to the directory where the output files will be saved")
args = p.parse_args()

SLIDE_FF_DIR = args.input
print(SLIDE_FF_DIR)
OUTPUT_DIR = args.output
OUTPUT_DIR.mkdir(parents=True, exist_ok=True) # create output dir if it doesn't exist


# %% looping through slides and creating objects

meta_file = [item for item in os.listdir(SLIDE_FF_DIR) if 'metadata_file' in item][0]
counts_file = [item for item in os.listdir(SLIDE_FF_DIR) if 'exprMat_file' in item][0]

# building the object
obj = sq.read.nanostring( # sq -> sc
    path=SLIDE_FF_DIR,
    counts_file=counts_file,
    meta_file=meta_file
)

# set cell_id as index
obj.obs = obj.obs.set_index('cell_id')

# modify the metadata file
obj.obs['slide_name'] = 'SLIDE0' + obj.obs['slide_ID'].astype(str)
obj.obs['fov_name'] = obj.obs['slide_name'].astype(str) + '_' + obj.obs['fov'].astype(str)

# %% saving metadata table
obj.obs.to_csv(
    OUTPUT_DIR / 'meta.tsv',
    sep='\t'
)

# saving the h5ad object
obj.write_h5ad(OUTPUT_DIR / 'adata.h5ad')

# free memory
del obj
gc.collect()
# %%
