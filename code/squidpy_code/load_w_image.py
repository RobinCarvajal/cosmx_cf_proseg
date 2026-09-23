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
    #path=sample_dir, #enable to add the images
    path=flat_file_dir,
    counts_file=counts_file, 
    meta_file=meta_file,
    fov_file=fov_file
)

end =  time.time()
print(f"Elapsed time: {end - start:.2f} seconds")

# %%
print(adata)


# %%

# This part of the analyis is similar to above
adata.var["Negative"] = adata.var_names.str.startswith("Negative")
adata.var["SystemControl"] = adata.var_names.str.startswith("SystemControl")
sc.pp.calculate_qc_metrics(adata, qc_vars=["Negative", "SystemControl"], inplace=True)
sc.pp.filter_cells(adata, min_counts=250)
sc.pp.filter_genes(adata, min_cells=1000)
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, inplace=True, exclude_highly_expressed=True)
sc.pp.log1p(adata)
sc.pp.pca(adata, n_comps=50)
sc.pp.neighbors(adata)
sc.tl.umap(adata, min_dist=0.2, spread=1)
sc.tl.leiden(adata)

# %%

adata_subset = adata[adata.obs.fov == "35"].copy()

# %%
sq.pl.spatial_segment(
    adata,
    color="leiden",
    seg_contourpx=20,
    seg_cell_id="cell_ID",
    library_key="fov",
    library_id="35",
    img=True,
    size=60,
    figsize = (4, 4),
    dpi = 200,
    #save = "fig-image1.png"
)


# %%

sq.gr.spatial_neighbors(adata_subset, coord_type="generic", delaunay=True)
sq.gr.spatial_autocorr(
    adata_subset,
    mode="moran",
    n_perms=100,
    n_jobs=1,
)
adata_subset.uns["moranI"].head(10)

sq.pl.spatial_segment(
    adata_subset,
    library_id="35",
    seg_cell_id="cell_ID",
    library_key="fov",
    color=["COX6C", "COL1A1", "KRT19"],
    size=60,
    img=False,
    figsize=(4, 4),
    dpi = 300,
    #save = "fig-nhood-Morans-I.png"
)


plt.show()


# %% save the object 

adata.write("F:/projects_work/cosmx_gray/data/slide_1.h5ad")