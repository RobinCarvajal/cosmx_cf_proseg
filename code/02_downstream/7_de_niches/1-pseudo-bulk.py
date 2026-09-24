# %% Libraries
import os
import scanpy as sc
#import scvi
import numpy as np
from pathlib import Path
import pandas as pd
from scipy.sparse import issparse

# %% Setting Paths
MAIN_DIR_NAME = "proseg_data"
MAIN_DIR = next(p for p in Path.cwd().parents if (p / MAIN_DIR_NAME).exists()) / MAIN_DIR_NAME
os.chdir(MAIN_DIR)

# %% Setting Seed
SEED_VALUE = 42
# set NumPy RNG for consistency
np.random.seed(SEED_VALUE)

# %% object versions
CUR_OBJ_V = 'annotated'
CUR_OBJ_PATH = MAIN_DIR / 'data' / 'comb-full' / 'h5ad' / f'comb-full-{CUR_OBJ_V}.h5ad'

# %% Load the unintegrated object
comb = sc.read_h5ad(CUR_OBJ_PATH)

# %% aggregate_counts function

def aggregate_counts(
    adata,
    celltype_col: str,
    samples_col: str,
    celltype_value: str,
    layer: str,
    shape: str = "samples-genes",  # "samples-genes" or "genes-samples"
):
    """
    Filter AnnData to one cell type and return a DataFrame of summed counts.

    Parameters
    ----------
    shape:
      - "samples-genes": rows = samples, columns = genes
      - "genes-samples": rows = genes, columns = samples
    """
    if shape not in {"samples-genes", "genes-samples"}:
        raise ValueError("shape must be one of: 'samples-genes', 'genes-samples'")

    # Filter to chosen cell type
    subset = adata[adata.obs[celltype_col] == celltype_value].copy()

    # Sample labels (one per cell)
    sample_ids = subset.obs[samples_col].astype(str).to_numpy()

    # Extract matrix from layer
    if layer not in subset.layers:
        raise ValueError(f"layer {layer!r} not found in adata.layers")

    X = subset.layers[layer]
    X = X.toarray() if issparse(X) else np.asarray(X)

    # cells × genes dataframe (index = sample id per cell)
    df_cells_genes = pd.DataFrame(X, index=sample_ids, columns=subset.var_names)

    # samples × genes
    df_samples_genes = df_cells_genes.groupby(level=0).sum()

    if shape == "samples-genes":
        return df_samples_genes
    else:  # "genes-samples"
        return df_samples_genes.T

# %% Exporting the pseudo-bulk tables per cell type

NICHES_VAR="novae_domains_11"

niche_list = comb.obs[NICHES_VAR].dropna().unique().tolist()

for niche in niche_list:
    print(f"Processing cell type: {niche}")

    table_dir = MAIN_DIR / 'results' / 'comb' / f'de_{NICHES_VAR}' / niche
    table_dir.mkdir(parents=True, exist_ok=True)

    pb_df = aggregate_counts(
        comb,
        celltype_col=NICHES_VAR,
        samples_col='sample_name',
        celltype_value=niche,
        layer='counts',
        shape='genes-samples',
        #shape='samples-genes'
    )

    # save pseudo-bulk dataframe
    out_dir = table_dir / f'pb_{niche.replace(" ", "_")}.csv'
    pb_df.to_csv(out_dir)

