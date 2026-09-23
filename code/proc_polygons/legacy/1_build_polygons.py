
# %% Libraries
import os
import scanpy as sc
import numpy as np
from pathlib import Path
import rapids_singlecell as rsc
import pandas as pd
import sys

# %% Setting Paths
MAIN_DIR_NAME = "cosmx_gray"
MAIN_DIR = next(p for p in Path.cwd().parents if (p / MAIN_DIR_NAME).exists()) / MAIN_DIR_NAME

# add it to sys.path to set it as root directory
sys.path.insert(0, str(MAIN_DIR))
# set as working directory
os.chdir(MAIN_DIR)

# %% custom libraries 
from utils.spatial_tools.helpers import map_annotations

# %% Setting Seed
SEED_VALUE = 42
# set NumPy RNG for consistency
np.random.seed(42)

# %% Load meta from comb object
meta = pd.read_csv('data/comb/meta/comb-meta.tsv', sep='\t', index_col=0)

# %% get all the cell names from the filtered comb object metadata
cell_names = meta.index.to_list()
len(cell_names)

# %% Load all polygon files from flat_files

FLAT_FILES_DIR = MAIN_DIR / 'data' / 'flat_files'

# list files in the directory 
flat_files = os.listdir(FLAT_FILES_DIR)
pattern = 'SLIDE'
# keep only files starting with 'SLIDE'
slide_dirs = [f for f in flat_files if f.startswith(pattern)]

filtered_poly_tbls = []

# loop 
for dir in slide_dirs:
    print(f"Processing directory: {dir}")
    dir_path = FLAT_FILES_DIR / dir
    # list all files in the directory
    all_files = os.listdir(dir_path)
    # filter for polygon files
    polygon_files = [f for f in all_files if f.endswith('-polygons.csv.gz')]
    
    for poly_file in polygon_files:
        print(f"  Processing file: {poly_file}")
        file_path = dir_path / poly_file
        # read the polygon file
        df = pd.read_csv(file_path)
        # filter for cells present in comb
        df_filtered = df[df['cell'].isin(cell_names)]
        print(f"    Original rows: {df.shape[0]}, Filtered rows: {df_filtered.shape[0]}")

        filtered_poly_tbls.append(df_filtered)

# %% Concatenate all filtered polygon tables
merged_poly_tbl = pd.concat(filtered_poly_tbls, ignore_index=True)
merged_poly_tbl.shape
# rename cell to cell_id
merged_poly_tbl = merged_poly_tbl.rename(columns={'cell': 'cell_id'})

# %% Save the merged and filtered polygon table
output_path = MAIN_DIR / 'data' / 'comb'/ 'polygons' 
output_path.mkdir(parents=True, exist_ok=True)

# %% Save the merged polygon table
merged_poly_tbl.to_csv(
    output_path / 'comb-polygons.csv.gz', 
    index=False, 
    compression='gzip'
)

# %%
