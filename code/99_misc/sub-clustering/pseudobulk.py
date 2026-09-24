# %% Libraries
from fontTools import subset
import os
import scanpy as sc
import numpy as np
from pathlib import Path
import rapids_singlecell as rsc
import pandas as pd
from scipy.sparse import issparse

# %% Setting Paths
MAIN_DIR_NAME = "cosmx_gray"
MAIN_DIR = next(p for p in Path.cwd().parents if (p / MAIN_DIR_NAME).exists()) / MAIN_DIR_NAME
os.chdir(MAIN_DIR)

# %% Setting Seed
SEED_VALUE = 42
# set NumPy RNG for consistency
np.random.seed(SEED_VALUE)

# %% object versions
CUR_OBJ_V = 'manual-niches'
CUR_OBJ_PATH = MAIN_DIR / 'data' / 'comb' / 'h5ad' / f'comb-{CUR_OBJ_V}.h5ad'

# %% Load the object
comb = sc.read_h5ad(CUR_OBJ_PATH)

# %% aggregate_counts function

import numpy as np
import pandas as pd
from scipy.sparse import issparse, csr_matrix


def aggregate_counts_multi(
    adata,
    celltype_col: str,
    samples_col: str,
    celltype_values: list,
    layer: str,
    shape: str = "genes-samples",  # "samples-genes" or "genes-samples"
    min_cells: int = 20,
    separator: str = "__",
    clean_names: bool = True,
):
    """
    Aggregate raw counts by sample × cell type.

    This is useful for pseudobulk DESeq2 analyses comparing:
      - secretory vs secretory-ciliated
      - ciliated vs secretory-ciliated
      - secretory vs ciliated

    Parameters
    ----------
    adata:
        AnnData object.

    celltype_col:
        Column in adata.obs containing cell type labels.

    samples_col:
        Column in adata.obs containing sample/donor IDs.

    celltype_values:
        Cell types to include, for example:
        ["secretory", "secretory-ciliated", "ciliated"]

    layer:
        Layer containing raw counts.

    shape:
        - "samples-genes": rows = pseudobulk samples, columns = genes
        - "genes-samples": rows = genes, columns = pseudobulk samples
          This is usually best for DESeq2 export.

    min_cells:
        Minimum number of cells required per sample × cell type pseudobulk.

    separator:
        Separator used in pseudobulk column names.

    clean_names:
        If True, replaces "-" with "_" in cell type names.
        This makes DESeq2 coefficient names easier to handle.

    Returns
    -------
    counts_df:
        Pseudobulk count matrix.

    meta_df:
        Metadata table for pseudobulk samples.
    """

    if shape not in {"samples-genes", "genes-samples"}:
        raise ValueError("shape must be one of: 'samples-genes', 'genes-samples'")

    if layer not in adata.layers:
        raise ValueError(f"Layer {layer!r} not found in adata.layers")

    # --------------------------------------------------
    # Filter to selected cell types
    # --------------------------------------------------
    mask = adata.obs[celltype_col].isin(celltype_values).to_numpy()

    obs = adata.obs.loc[mask, [samples_col, celltype_col]].copy()

    X = adata.layers[layer][mask, :]

    if issparse(X):
        X = X.tocsr()
    else:
        X = np.asarray(X)

    # --------------------------------------------------
    # Clean labels
    # --------------------------------------------------
    obs[samples_col] = obs[samples_col].astype(str)
    obs[celltype_col] = obs[celltype_col].astype(str)

    if clean_names:
        obs[celltype_col] = (
            obs[celltype_col]
            .str.replace("-", "_", regex=False)
            .str.replace(" ", "_", regex=False)
        )

    # Create sample × cell type pseudobulk IDs
    obs["pseudobulk_id"] = (
        obs[samples_col] + separator + obs[celltype_col]
    )

    # --------------------------------------------------
    # Aggregate counts efficiently
    # --------------------------------------------------
    group_codes, group_names = pd.factorize(obs["pseudobulk_id"], sort=True)

    n_groups = len(group_names)
    n_cells = obs.shape[0]

    design = csr_matrix(
        (
            np.ones(n_cells, dtype=np.float32),
            (group_codes, np.arange(n_cells))
        ),
        shape=(n_groups, n_cells)
    )

    pb_counts = design @ X

    if issparse(pb_counts):
        pb_counts = pb_counts.toarray()

    counts_samples_genes = pd.DataFrame(
        pb_counts,
        index=group_names,
        columns=adata.var_names
    )

    # --------------------------------------------------
    # Metadata
    # --------------------------------------------------
    meta_df = (
        obs
        .drop_duplicates("pseudobulk_id")
        .set_index("pseudobulk_id")
        .loc[group_names]
        .copy()
    )

    meta_df = meta_df.rename(
        columns={
            samples_col: "sample_name",
            celltype_col: "cell_type"
        }
    )

    meta_df["n_cells"] = np.bincount(group_codes)

    # --------------------------------------------------
    # Remove pseudobulks with too few cells
    # --------------------------------------------------
    keep = meta_df["n_cells"] >= min_cells

    counts_samples_genes = counts_samples_genes.loc[keep]
    meta_df = meta_df.loc[keep]

    # --------------------------------------------------
    # Return desired orientation
    # --------------------------------------------------
    if shape == "samples-genes":
        counts_df = counts_samples_genes
    else:
        counts_df = counts_samples_genes.T

    return counts_df, meta_df

# %% Make a subset only for secretory, ciliated and secretory-ciliated cells

celltypes_to_keep = [
    "secretory",
    "secretory-ciliated",
    "ciliated"
]

pb_counts, pb_meta = aggregate_counts_multi(
    adata=comb,
    celltype_col="ann_lvl_3",      # change if yours is different
    samples_col="sample_name",
    celltype_values=celltypes_to_keep,
    layer="counts",               # use your raw counts layer
    shape="genes-samples",
    min_cells=20
)

# %% test

pb_counts.head()
