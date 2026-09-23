import plotnine as p9
import pandas as pd


def dim_plot(
        adata,
        reduction, 
        group_by, 
        split_by=None, 
        pt_size=0.1, 
        add_labels=True, 
        fig_size=(12, 12), 
        equal_coords=True,
        legend_name=None,
        legend_position='right',
        legend_font_size=12,
        rasterize=True,
        ):
    """
    Make a UMAP plot with optional faceting and Seurat-style group labels.
    
    Parameters
    ----------
    adata : AnnData
        AnnData object containing single-cell data.
    reduction : str
        The embedding key to use from `adata.obsm`, e.g. "X_umap" or "X_pca".
    group_by : str
        Column name in `adata.obs` to color the points by (e.g. clusters or cell types).
    split_by : str, optional
        Column name in `adata.obs` to facet the plot (create multiple panels).
    pt_size : float
        Size of the points in the scatter plot.
    add_labels : bool
        Whether to add group labels (cluster-style labels) at the median position.
    fig_size : tuple
        Figure size in inches (width, height).
    """

    # Extract the embedding (e.g. UMAP coordinates) from `adata.obsm`
    df = adata.obsm[reduction]

    # Define column names for plotting
    x_col = f'{reduction}_1'
    y_col = f'{reduction}_2'

    # Create a DataFrame with x and y coordinates
    df = pd.DataFrame(df, columns=[x_col, y_col])

    # If facetting by a variable, add that metadata column to the DataFrame
    if split_by:
        df[split_by] = adata.obs[split_by].values

    # Always add the main grouping column for coloring
    if group_by:
        df[group_by] = adata.obs[group_by].values

    # Create the base ggplot scatter plot
    p = (
        p9.ggplot(df, p9.aes(x=x_col, y=y_col, color=group_by))
        + p9.geom_point(shape='.', size=pt_size, raster=rasterize)
        + p9.theme(
            figure_size=fig_size,
            legend_position=legend_position,
            legend_text=p9.element_text(size=legend_font_size),
            legend_title=p9.element_text(size=legend_font_size),
            legend_key_size=legend_font_size,
        )
        + p9.guides(color=p9.guide_legend(override_aes={'size': legend_font_size}))
    )

    # If facetting is requested, split plot into sub-panels by `split_by`
    if split_by:
        p += p9.facet_wrap(f'~{split_by}')

    # If adding labels, compute centroids and overlay text labels
    if add_labels:
        centroids = (
            df.groupby(group_by)[[x_col, y_col]]
            .median()             # take median x,y of each group
            .reset_index()        # turn group labels into a column
        )

        # Add labels at the group centroids
        p += p9.geom_label(
            centroids,
            p9.aes(x=x_col, y=y_col, label=group_by, fill=group_by),
            color="black",      # label text color
            size=legend_font_size,             # label font size
            ha="center",         # center align text horizontally
            show_legend=False,  # don't show legend for labels
            fontweight="bold",       # make label text bold
        )

    # Set equal aspect ratio
    if equal_coords:
        p += p9.coord_equal()

    # legend title
    if legend_name:
        p += p9.labs(color=legend_name)

    # Return the ggplot object
    return p
