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

# output directory 
OUT_DIR = MAIN_DIR / 'results' / 'comb-full' / 'de_ct'

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

LABELS_KEY = 'ann_lvl_3_refined'
celltypes = [
    celltype
    for celltype in comb.obs[LABELS_KEY].unique().tolist()
]

for ct in celltypes:
    print(f"\nProcessing cell type: {ct}")

    try:
        pb_df = aggregate_counts(
            comb,
            celltype_col=LABELS_KEY,
            samples_col="donor",
            celltype_value=ct,
            layer="counts",
            shape="genes-samples",
        )

        # Skip empty pseudobulk tables
        if pb_df.empty or pb_df.shape[1] == 0:
            print(f"Skipping {ct}: pseudobulk table is empty")
            continue

        counts_array = pb_df.to_numpy()

        if not np.isfinite(counts_array).all():
            print(f"Skipping {ct}: contains NaN or infinite values")
            continue

        if (counts_array < 0).any():
            print(f"Skipping {ct}: contains negative counts")
            continue

        if not np.allclose(counts_array, np.round(counts_array)):
            print(f"Skipping {ct}: contains non-integer counts")
            continue

        # Convert to integer raw counts
        pb_df = pb_df.astype(np.int64)

        # Keep genes with mean raw count >= 1
        pb_df = pb_df.loc[pb_df.mean(axis=1) >= 1]

        # Skip if filtering removed every gene
        if pb_df.empty:
            print(f"Skipping {ct}: no genes remain after filtering")
            continue

        table_dir = OUT_DIR / ct
        table_dir.mkdir(parents=True, exist_ok=True)

        out_path = table_dir / f'pb_{ct.replace(" ", "_")}.csv'
        pb_df.to_csv(out_path)

        print(
            f"Saved {ct}: "
            f"{pb_df.shape[0]} genes × {pb_df.shape[1]} samples"
        )

    except Exception as error:
        print(f"Skipping {ct} because of error: {error}")
        continue

