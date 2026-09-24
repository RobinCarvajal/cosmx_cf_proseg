
# %% Libraries

import scanpy as sc
import pandas as pd
from qc_utils import run_qc, qc_filter
from pathlib import Path
import os
import gc

# %% setting paths
MAIN_DIR_NAME = "cosmx_gray"
MAIN_DIR = next(p for p in Path.cwd().parents if (p / MAIN_DIR_NAME).exists()) / MAIN_DIR_NAME
os.chdir(MAIN_DIR)

# %% extract slide names from the slides directory

SAMPLES_DIR = Path('data/samples')

sample_names = [
    d.name for d in SAMPLES_DIR.iterdir() 
    if d.is_dir() and (
        "slide" in d.name.lower() or
        "ctrl" in d.name.lower() or
        "pwcf" in d.name.lower()
    )
]
print(f"Samples found: {sample_names}")

OBJ_VERSION = 'v0'
NEW_OBJ_VERSION = 'v1'

# %% loop through samples and run qc

for sample in sample_names:
    
    print(f'Processing {sample}')

    # load the sample object
    sample_path = SAMPLES_DIR / f'{sample}/h5ad/{sample}_{OBJ_VERSION}.h5ad'
    sample_obj = sc.read_h5ad(sample_path)

    # run qc
    sample_obj = run_qc(sample_obj)
    print(sample_obj)

    # Create a quick summary table of all QC outcomes (counts of pass/fail)
    qc_cols = [
        'qc_cell_counts', 
        'qc_fov_counts', 
        'qc_complexity', 
        'qc_negprobes', 
        'qc_borders', 
        'qc_pass'
    ]
    
    summary = sample_obj.obs[qc_cols].apply(pd.Series.value_counts)
    # save the summary table
    summary_path = SAMPLES_DIR / f'{sample}/tables/{sample}_qc_summary.tsv'
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(summary_path, sep='\t')
    print(f'QC summary saved at {summary_path}')

    # filter the object
    sample_obj = qc_filter(sample_obj)
    print(sample_obj)

    # save the filtered object
    sample_path_filt = SAMPLES_DIR / f'{sample}/h5ad/{sample}_{NEW_OBJ_VERSION}.h5ad'
    sample_obj.write_h5ad(sample_path_filt)

    print(f'Object saved at {sample_path_filt}')

    # free memory
    del sample_obj
    gc.collect()

# %%
