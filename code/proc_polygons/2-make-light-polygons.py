# %%
import pandas as pd 
import os
from pathlib import Path

# %% Set Main Directory
MAIN_DIR_NAME = "cosmx_gray"
MAIN_DIR = next(p for p in Path.cwd().parents if (p / MAIN_DIR_NAME).exists()) / MAIN_DIR_NAME

# add it to sys.path to set it as root directory
os.chdir(MAIN_DIR)

# %% load meta data
polygons = pd.read_parquet('data/comb/polygons/comb-polygons.parquet')

# %% save a lighter version only with cell_id and global coords
polygons_light = polygons[["cell_id","x_global_px","y_global_px"]]
polygons_light = polygons_light.rename(columns={"x_global_px": "x", "y_global_px": "y"})

# %% save as parquet
polygons_light.to_parquet("data/comb/polygons/comb-polygons-light.parquet", index=False)