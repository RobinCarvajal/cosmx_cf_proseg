import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import scanpy as sc
from anndata import AnnData
from scipy.sparse import issparse

def silverman_bandwidth(data):
    """
    Compute bandwidth using Silverman's rule of thumb:
    h = 0.9 * min(std, IQR/1.34) * n^(-1/5)
    """
    n = len(data)
    std = np.std(data, ddof=1)
    iqr = np.subtract(*np.percentile(data, [75, 25]))
    return 0.9 * min(std, iqr / 1.34) * n ** (-1/5)

def wkde2d(x, y, w, adjust=1, n=100, lims=None):
    """
    Weighted 2D Kernel Density Estimation (KDE)

    Args:
      x, y: 1D arrays of coordinates for data points.
      w:    1D array of weights (e.g., gene expression values).
      adjust: Bandwidth scaling factor.
      n:      Grid size (number of points along x and y).
      lims:   [xmin, xmax, ymin, ymax]. If None, computed from data.

    Returns:
      gx, gy: 1D arrays of grid coordinates along x and y.
      z:      2D density matrix of shape (n, n).
    """
    x = np.asarray(x)
    y = np.asarray(y)
    w = np.asarray(w)
    if len(x) != len(y) or len(x) != len(w):
        raise ValueError("x, y, and w must have the same length")
    if lims is None:
        lims = [np.min(x), np.max(x), np.min(y), np.max(y)]
    
    h_x = silverman_bandwidth(x) * adjust
    h_y = silverman_bandwidth(y) * adjust

    gx = np.linspace(lims[0], lims[1], n)
    gy = np.linspace(lims[2], lims[3], n)

    # Standardized differences between grid points and data points
    ax = (gx[:, None] - x[None, :]) / h_x  # shape: (n, N)
    ay = (gy[:, None] - y[None, :]) / h_y  # shape: (n, N)

    pdf_ax = norm.pdf(ax)
    pdf_ay = norm.pdf(ay)

    # Broadcast weights across rows
    w_mat = np.tile(w, (n, 1))

    # Apply kernel values
    A = pdf_ax * w_mat
    B = pdf_ay * w_mat

    # Density for each grid point (gx, gy), aggregating weighted contributions
    z = np.dot(A, B.T) / (np.sum(w) * h_x * h_y)
    return gx, gy, z

def nebulosa_density(adata, coord_key, gene, adjust=1, n=100, lims=None, cmap='viridis', show=False, pt_size=1):
    """
    For an AnnData object, compute weighted 2D KDE using a 2D embedding
    stored in .obsm and a gene's expression values, then (optionally) plot.

    Args:
      adata:     AnnData object.
      coord_key: Key in adata.obsm for 2D coordinates (e.g., "X_umap"), shape (n_cells, 2).
      gene:      Gene name used as weights (must exist in adata.var_names).
      adjust:    Bandwidth scaling factor (default 1).
      n:         Grid size (default 100).
      lims:      [xmin, xmax, ymin, ymax]. If None, computed from data.
      cmap:      Matplotlib colormap for plotting (default 'viridis').
      show:      If True, call plt.show() to display the figure.

    Returns:
      If show is False:
        densities: 1D array of KDE values mapped to each original cell.
      If show is True:
        fig, ax: Matplotlib Figure and Axes objects.
    """
    # Extract 2D coordinates
    if coord_key not in adata.obsm.keys():
        raise KeyError(f"{coord_key} is not present in adata.obsm")
    coords = adata.obsm[coord_key]
    if coords.shape[1] < 2:
        raise ValueError("The selected coordinates must have at least two columns")
    coords = coords[:, :2]
    x = coords[:, 0]
    y = coords[:, 1]
    
    # Check that the gene exists
    if gene not in adata.var_names:
        raise KeyError(f"Gene {gene} is not present in adata.var_names")
    # Get expression values as a 1D array
    expr = adata[:, gene].X
    if issparse(expr):
        expr = expr.toarray().flatten()
    else:
        expr = np.array(expr).flatten()
    
    # Compute weighted KDE
    gx, gy, z = wkde2d(x, y, expr, adjust=adjust, n=n, lims=lims)
    
    # Map each original point onto the grid to obtain density values
    ix = np.digitize(x, gx) - 1
    iy = np.digitize(y, gy) - 1
    ix = np.clip(ix, 0, len(gx)-1)
    iy = np.clip(iy, 0, len(gy)-1)
    densities = z[ix, iy]
    # if not show:
    #     return densities
    
    # Scatter plot colored by density
    plt.figure(figsize=(8, 6))
    sc = plt.scatter(x, y, c=densities, cmap=cmap, s=pt_size)
    plt.colorbar(sc, label='Weighted KDE Density')
    plt.xlabel('Dimension 1')
    plt.ylabel('Dimension 2')
    plt.title(f'Weighted KDE Density for {gene}')

    # 🔑 Fix scaling so one unit on x = one unit on y
    plt.gca().set_aspect('equal', adjustable='box')

    fig = plt.gcf()
    if not show:
        plt.close(fig)

    return fig


def nebulosa_joint_density(
    adata, coord_key, genes, 
    adjust=1, n=100, lims=None, cmap='viridis', 
    show=False, pt_size=1, method="product", threshold=0
):
    """
    For an AnnData object, compute weighted 2D KDE using a 2D embedding
    stored in .obsm and the expression values of multiple genes.

    Args:
      adata:     AnnData object.
      coord_key: Key in adata.obsm for 2D coordinates (e.g., "X_umap").
      genes:     List of gene names (must exist in adata.var_names).
      adjust:    Bandwidth scaling factor.
      n:         Grid size.
      lims:      [xmin, xmax, ymin, ymax]. If None, computed from data.
      cmap:      Matplotlib colormap for plotting.
      show:      If True, display the figure.
      pt_size:   Point size in scatter plot.
      method:    How to combine expressions across genes.
                 Options: "product", "sum", "mean", "boolean".
      threshold: Expression cutoff for "boolean" method (default 0).

    Returns:
      fig: Matplotlib Figure.
    """
    if coord_key not in adata.obsm:
        raise KeyError(f"{coord_key} not found in adata.obsm")

    coords = adata.obsm[coord_key][:, :2]
    x, y = coords[:, 0], coords[:, 1]

    # ensure genes exist
    missing = [g for g in genes if g not in adata.var_names]
    if missing:
        raise ValueError(f"Genes not found in adata.var_names: {missing}")

    # extract expression matrix for selected genes
    expr_mat = adata[:, genes].X.toarray()

    # combine according to method
    if method == "product":
        expr = np.prod(expr_mat, axis=1)
    elif method == "sum":
        expr = np.sum(expr_mat, axis=1)
    elif method == "mean":
        expr = np.mean(expr_mat, axis=1)
    elif method == "boolean":
        expr = np.all(expr_mat > threshold, axis=1).astype(int)
    else:
        raise ValueError("method must be one of: 'product', 'sum', 'mean', 'boolean'")

    # Compute weighted KDE
    gx, gy, z = wkde2d(x, y, expr, adjust=adjust, n=n, lims=lims)

    # Map each cell back to grid
    ix = np.clip(np.digitize(x, gx) - 1, 0, len(gx)-1)
    iy = np.clip(np.digitize(y, gy) - 1, 0, len(gy)-1)
    densities = z[ix, iy]

    # Plot
    plt.figure(figsize=(8, 6))
    sc = plt.scatter(x, y, c=densities, cmap=cmap, s=pt_size)
    plt.colorbar(sc, label='Weighted KDE Density')
    plt.xlabel('Dimension 1')
    plt.ylabel('Dimension 2')
    plt.title(f'Weighted KDE Density for {", ".join(genes)} [{method}]')

    # 🔑 Fix scaling so one unit on x = one unit on y
    plt.gca().set_aspect('equal', adjustable='box')

    fig = plt.gcf()
    if not show:
        plt.close(fig)

    return fig



