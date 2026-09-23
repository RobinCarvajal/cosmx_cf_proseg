# %% load libraries
from pathlib import Path
import os
import scanpy as sc
import squidpy as sq
import gc

# %% setting paths
MAIN_DIR_NAME = "cosmx_gray"
MAIN_DIR = next(p for p in Path.cwd().parents if (p / MAIN_DIR_NAME).exists()) / MAIN_DIR_NAME
os.chdir(MAIN_DIR)

# %% loading flat files

FLATFILES_DIR = Path('data/flat_files')

# dirs inside flatfiles files
slide_names = [item for item in os.listdir(FLATFILES_DIR) if os.path.isdir(FLATFILES_DIR / item)]

OBJ_VERSION = 'v0'

# %% looping through slides and creating objects

for slide in slide_names:
    
    print(slide)

    # creating output dirs
    h5ad_dir = Path( f'data/slides/{slide}/h5ad')
    tbls_dir = Path( f'data/slides/{slide}/tables')
    h5ad_dir.mkdir(parents=True, exist_ok=True)
    tbls_dir.mkdir(parents=True, exist_ok=True)

    slide_ff_dir = FLATFILES_DIR / slide
    meta_file = [item for item in os.listdir(slide_ff_dir) if 'metadata_file' in item][0]
    counts_file = [item for item in os.listdir(slide_ff_dir) if 'exprMat_file' in item][0]

    # building the object
    obj = sq.read.nanostring(
        path=slide_ff_dir,
        counts_file=counts_file,
        meta_file=meta_file
    )

    # set cell_id as index
    obj.obs = obj.obs.set_index('cell_id')

    # modify the metadata file
    obj.obs['slide_name'] = 'SLIDE0' + obj.obs['slide_ID'].astype(str)
    obj.obs['fov_name'] = obj.obs['slide_name'].astype(str) + '_' + obj.obs['fov'].astype(str)

    # saving metadata table
    obj.obs.to_csv(
        tbls_dir / f'{slide}_meta.tsv',
        sep='\t'
    )

    # saving the h5ad object
    obj.write_h5ad(h5ad_dir / f'{slide}_{OBJ_VERSION}.h5ad')

    # free memory
    del obj
    gc.collect()
# %%
