
# %% Load libs 

import scanpy as sc
import pandas as pd
import numpy as np
import os 
import anndata as ad
from pathlib import Path

# %% Load h5ad objects

MAIN_DIR_NAME = "cosmx_gray"
MAIN_DIR = next(p for p in Path.cwd().parents if (p / MAIN_DIR_NAME).exists()) / MAIN_DIR_NAME
os.chdir(MAIN_DIR)

# Define the samples directory
samples_dir = "data/samples"

# get the sample names
sample_names_lst = os.listdir(samples_dir)
# keep only names that have "sample" in them
sample_names_lst = [name for name in sample_names_lst if "sample" in name]

# list to store sample objects
sample_obj_dict = {}

# loop 
for sample in sample_names_lst:

    # get the h5ad path
    h5ad_path = os.path.join(samples_dir, sample, "h5ad", f"{sample}_v1.h5ad")

    # read the h5ad file
    sample_obj = sc.read(h5ad_path)

    # store in list
    sample_obj_dict[sample] = sample_obj

# %% Define high quality cell (hq_cell) based on QC metrics

for sample in sample_names_lst:

    # define the output path
    sample_obj = sample_obj_dict[sample]

    # define new column
    sample_obj.obs["hq_cell"] = np.where(
        (sample_obj.obs["qc_complexity"] == "pass") &
        (sample_obj.obs["qc_prop_negprobes"] == "pass") &
        (sample_obj.obs["qc_count"] == "pass") &
        (sample_obj.obs["qc_border_cell"] == 0) & # 0 represents False
        (sample_obj.obs["qc_solidity"] == "pass") &
        (sample_obj.obs["qc_fov"] == "pass"),
        True,
        False
    )
    print(sample)
    print(sample_obj.obs["hq_cell"].value_counts())


# %% Keep only high quality cells

for sample in sample_names_lst:

    # define the output path
    sample_obj = sample_obj_dict[sample]

    # keep only high quality cells
    sample_obj = sample_obj[sample_obj.obs["hq_cell"] == True]

    # save the object back to the list 
    sample_obj_dict[sample] = sample_obj
    
    print(sample)
    print(sample_obj.obs["hq_cell"].value_counts())

# %% Merge all samples into one object 

comb = ad.concat(sample_obj_dict)
comb.obs_names_make_unique()
print(comb.obs["sample_name"].value_counts())
comb

# %%
sc.pl.violin(
    comb,
    ["nCount_RNA", "nFeature_RNA"],
    jitter=0.1,
    size=0.09,
    multi_panel=True,
)
# %%

sc.pl.scatter(
    comb, 
    "nCount_RNA", 
    "nFeature_RNA", 
    color="qc_area",
    size = 3
)

# %%

sc.pp.filter_genes(comb, min_cells=3)

# %% create folder strcutre and save the object

from pathlib import Path

comb_dir = Path('data/comb')
subdirs = ['h5ad', 'meta']

# create paths
paths = [comb_dir / sd for sd in subdirs]

# make directories (including comb itself)
for p in paths:
    p.mkdir(parents=True, exist_ok=True)

print(paths)
# [PosixPath('comb/h5ad'), PosixPath('comb/meta')]


# %%
h5ad_path = "data/comb/h5ad/comb_v0.h5ad"
comb.write(h5ad_path)

# %% save meta as tsv
meta = comb.obs.copy()
# save 
META_PATH = MAIN_DIR / 'data/comb/meta' / 'comb-meta.tsv'
meta.to_csv(META_PATH, sep='\t', index=True, header=True)
# %%
