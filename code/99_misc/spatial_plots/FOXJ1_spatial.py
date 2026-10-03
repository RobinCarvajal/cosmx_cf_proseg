# %% Libraries 
import os
import sys
import pandas as pd
import plotnine as p9
from pathlib import Path
import pypalettes as pp
import matplotlib.colors as mcolors
import scanpy as sc
import numpy as np

# %% Set Main Directory
MAIN_DIR_NAME = "cosmx_cf_proseg"
MAIN_DIR = next(p for p in Path.cwd().parents if (p / MAIN_DIR_NAME).exists()) / MAIN_DIR_NAME

# add it to sys.path to set it as root directory
os.chdir(MAIN_DIR)

# %% Results Directories
PLOTS_DIR = Path("results/comb/spatial_polygons/FOXJ1")
print(PLOTS_DIR)

# %% CRATE code

# imports
import os
import sys
import pandas as pd
import plotnine as p9
from pathlib import Path
import pypalettes as pp
import matplotlib.colors as mcolors
import scanpy as sc

# create crate
def create_crate(adata_path, polygons_path, id_col):
    
    adata = sc.read_h5ad(adata_path)
    polygons = pd.read_parquet(polygons_path)

    crate = {
        "adata": adata,
        "polygons": polygons,
        "id_col": id_col,
    }

    return crate

# subset crate
def subset_crate(crate, filter_mask):

    id_col = crate["id_col"]
    adata = crate["adata"]
    polygons = crate["polygons"]

    cells = adata.obs.index[filter_mask]
    sub_adata = adata[cells, :].copy()
    sub_polygons = polygons.loc[polygons[id_col].isin(cells)].copy()

    sub_crate = {
        "adata": sub_adata,
        "polygons": sub_polygons,
        "id_col": id_col
    }
    
    return sub_crate

# plot_polygons
def plot_polygons(
    crate, 
    ann_var,
    fig_size: tuple = (20, 20),
    colors=None,
    layer="counts",
    cmap="viridis",
    background_color:str="black"
):
    import numpy as np
    import pandas as pd
    import plotnine as p9
    from pandas.api.types import is_numeric_dtype

    # if ann_var not in crate["adata"].obs.columns look into adata.var_names
    if ann_var not in crate["adata"].obs.columns:
        if ann_var in crate["adata"].var_names:
            # if it's in var_names, we need to map it to obs using the layers
            expr_values = pd.DataFrame(
                crate["adata"][:, ann_var].layers[layer].toarray(),
                index=crate["adata"].obs_names,
                columns=[ann_var]
            )

            # log1p transform
            expr_values = np.log1p(expr_values)

            # safer than merge into whole obs
            crate["adata"].obs[ann_var] = expr_values.loc[crate["adata"].obs_names, ann_var].values
        else:
            raise ValueError(
                f"Annotation variable {ann_var} not found in adata.obs or adata.var_names."
            )

    id_col = crate["id_col"]
    meta_df = crate["adata"].obs[[ann_var]].copy()
    poly_df = crate["polygons"]

    #### map ann_var to polygons ####
    ann_poly = pd.merge(poly_df, meta_df, on=id_col, how="left")

    #### choose scale ####
    if is_numeric_dtype(ann_poly[ann_var]):
        fill_scale = (
            p9.scale_fill_cmap(name=cmap)
            if colors is None
            else p9.scale_fill_gradientn(colors=colors)
        )
    else:
        fill_scale = p9.scale_fill_manual(values=colors)

    #### plot ####
    p = (
        p9.ggplot(ann_poly)
        + p9.geom_map(
            p9.aes(
                geometry="geometry",
                fill=ann_var,
            ),
            color="white",
            size=0.1,
        )
        + p9.coord_fixed(1)
        + fill_scale
        + p9.theme(
            axis_line=p9.element_blank(),
            axis_text=p9.element_blank(),
            axis_ticks=p9.element_blank(),
            axis_title=p9.element_blank(),
            panel_background=p9.element_rect(fill=background_color),
            panel_grid_major=p9.element_blank(),
            panel_grid_minor=p9.element_blank(),
            figure_size=fig_size,
        )
    )

    return p

# %% Load adata and polygons

# Paths
OBJ_V = "cellcharter"
adata_path = Path(f"data/comb/h5ad/comb-{OBJ_V}.h5ad")
polygons_path = Path("data/comb/polygons/comb-polygons-light.parquet")
id_col = "cell_id"

# create crate
crate = create_crate(adata_path, polygons_path, id_col)

# %% get sample names
sample_names = crate["adata"].obs["sample_name"].unique()
# %% Set the annotation variable
VARS = ["FOXJ1"]

# %% Palette

# # Load the colormap
# cmap = pp.load_cmap("alphabet")
# # Get N discrete colors as hex
# # use leiden_scVI because it more clusters, keeps colouring consistent
# n = len(crate["adata"].obs[ANN_VAR].unique())
# palette = [cmap(i / (n-1)) for i in range(n)]
# # Convert RGBA -> hex
# palette = [mcolors.to_hex(c) for c in palette]


# %% Plot polygons

for sample in sample_names:

    print(sample)

    for ANN_VAR in VARS:

        print(ANN_VAR)

        # filter mask 
        filter_mask = crate["adata"].obs["sample_name"] == sample

        # subset crate
        sub_crate = subset_crate(crate, filter_mask)

        p = plot_polygons(sub_crate, ann_var=ANN_VAR, fig_size=(20,20))

        for ext in ['png', 'svg']:
            p.save(
                    filename= PLOTS_DIR / f'{sample}_{ANN_VAR}.{ext}',
                    dpi=900,
                    units='in'
                )