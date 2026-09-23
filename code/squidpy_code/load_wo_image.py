# -*- coding: utf-8 -*-
"""
Created on Mon May 26 15:01:52 2025

@author: robin
"""
# %% libraries

from pathlib import Path
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import squidpy as sq

# %% creating seurat object

import time 
start = time.time()

sample_dir = "/Volumes/robin_work/cosmx_gray/raw_files/Raw_Files_Slide_1/sample_dir_formatted"
flat_file_dir = "/Volumes/robin_work/cosmx_gray/raw_files/Raw_Files_Slide_1/flatFiles/NCLRGray08052025"
meta_file = [item for item in os.listdir(sample_dir) if 'metadata_file' in item][0]
counts_file = [item for item in os.listdir(sample_dir) if 'exprMat_file' in item][0]
fov_file = [item for item in os.listdir(sample_dir) if 'fov_positions_file' in item][0]
 
fov_df = pd.read_csv(os.path.join(sample_dir, fov_file))
if 'FOV' in fov_df.columns:
  print("Refactoring file to older format.")
  # Rename 'FOV' column to 'fov'
  fov_df.rename(columns={'FOV': 'fov'}, inplace=True)
  # have fov_file reference the new, formatted file and write it
  fov_file = os.path.join(sample_dir,'fov_positions_formatted.csv')
  fov_df.to_csv(fov_file, index=False)

adata = sq.read.nanostring(
    path=flat_file_dir,
    counts_file=counts_file, 
    meta_file=meta_file,
    fov_file=fov_file
)

end =  time.time()
print(f"Elapsed time: {end - start:.2f} seconds")

# %%
print(adata)

# %% can get a plot with this code usign only centroids

sq.pl.spatial_scatter(adata, shape='hex',color='Max.PanCK', library_key='fov', img=False, library_id=['1','2'], spatial_key='spatial', size=100)






















