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
os.chdir(MAIN_DIR)

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
meta.reset_index(inplace=True)  # reset index to have cell_id as a column

# %% Load the polygons file 
transcripts = pd.read_csv('data/comb/transcripts/comb-transcripts.csv')

# %% map the annotations in the meta to the transcripts file
transcripts_ann = map_annotations(
    vertices_df=transcripts,
    meta_df=meta,  # reset index to have cell_id as a column
    map_by='cell_id',
    cols_to_map=['sample_name']
)

# %% Save the merged and filtered transcript table
output_path = MAIN_DIR / 'data' / 'comb'/ 'transcripts' 
output_path.mkdir(parents=True, exist_ok=True)

# %% Save the merged transcript table
transcripts_ann.to_csv(
    output_path / 'comb-transcripts-annotated.csv', 
    index=False,
)
# %% Save the merged transcript table (compressed)
transcripts_ann.to_csv(
    output_path / 'comb-transcripts-annotated.csv.gz', 
    index=False, 
    compression='gzip'
)

# %%
