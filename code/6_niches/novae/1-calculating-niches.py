# %%
# %% Libraries
import os
import scanpy as sc
import numpy as np
from pathlib import Path
import novae
import time

# %%
# %% Setting Paths
MAIN_DIR_NAME = "cosmx_gray"
MAIN_DIR = next(p for p in Path.cwd().parents if (p / MAIN_DIR_NAME).exists()) / MAIN_DIR_NAME
os.chdir(MAIN_DIR)

# %% Setting Seed
SEED_VALUE = 42
# set NumPy RNG for consistency
np.random.seed(SEED_VALUE)

# %% object versions
CUR_OBJ_V = 'annotated'
CUR_OBJ_PATH = MAIN_DIR / 'data' / 'comb-full' / 'h5ad' / f'comb-full-{CUR_OBJ_V}.h5ad'
NEW_OBJ_V = 'annotated'
NEW_OBJ_PATH = MAIN_DIR / 'data' / 'comb-full' / 'h5ad' / f'comb-full-{NEW_OBJ_V}.h5ad'

# %% Load the unintegrated object
# use adata isntead of comb for the object only in this script
adata = sc.read_h5ad(CUR_OBJ_PATH)

# %% Transform spatial coords to um

# Need to transform to um from pixels. Information on measure transformation here:
# https://nanostring.com/wp-content/uploads/2023/09/SMI-ReadMe-BETA_humanBrainRelease.html

# save orginal spatial coordinates (i think they are local coordinates in px per fov)
adata.obsm["spatial_orig"] = adata.obsm["spatial"].copy()
# convert spatial from px to microns
adata.obsm["spatial_px"] = adata.obsm["X_spatial"].copy()
adata.obsm["spatial_um"] = adata.obsm["spatial_px"] * 0.12  # 0.12 to convert to um
adata.obsm["spatial"] = adata.obsm["spatial_um"]

# %% Calc spatial neighbours

# if you have multiple samples in the same adata object, specify `slide_key`
novae.spatial_neighbors(adata, radius=100, slide_key='sample_name', coord_type='generic') # this for multiple slides
# here used 'sample_name' because there is multiple samples per slide

# %%
novae.plot.connectivities(adata)

# %% [markdown]
# Use pretrained model

# %%
model = novae.Novae.from_pretrained("MICS-Lab/novae-human-0")
model

# %% Train model 
start_time = time.time()  # Record the start time

model.fine_tune(
    adata,
    accelerator='cuda', 
    num_workers=8,
    max_epochs=50,
) 

end_time = time.time()    # Record the end time
elapsed_time = end_time - start_time
print(f"Execution time: {elapsed_time:.4f} seconds")

# %% Compute representations

start_time = time.time()

model.compute_representations(
    adata, 
    zero_shot=False, 
    accelerator='cuda', 
    num_workers=8
)

end_time = time.time()
elapsed_time = end_time - start_time
print(f"Execution time: {elapsed_time:.4f} seconds")

# %% save in the metadata the spatial domain labels up to 12 levels
for i in range(1, 12):
  model.assign_domains(adata, i)

# %% Save the model
novae_model_path = MAIN_DIR / 'data' / 'comb-full' / 'novae_model'
# Create directory
novae_model_path.mkdir(parents=True, exist_ok=True)
# save
model.save_pretrained(save_directory=novae_model_path)

# %% Save the object with spatial domains annotations
adata.write_h5ad(NEW_OBJ_PATH)


