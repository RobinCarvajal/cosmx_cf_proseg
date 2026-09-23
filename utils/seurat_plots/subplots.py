import math
import numpy as np
import matplotlib.pyplot as plt

def _safe_draw(gg):
    """
    Draw a plotnine ggplot and return its matplotlib Figure.
    If patchworklib monkey-patched plotnine, reload to restore draw().
    """
    try:
        return gg.draw()
    except TypeError:
        # Recover from patchworklib's modified draw()
        import importlib, plotnine.ggplot as _ggmod
        importlib.reload(_ggmod)
        return gg.draw()

def render_ggplot_to_ax(gg, ax, dpi=200, close=True):
    """
    Render a plotnine ggplot into an existing Matplotlib Axes using imshow.
    - Forces a canvas draw for crisp text.
    - Copies pixels, then closes the temporary Figure to avoid leaks.
    """
    fig_tmp = _safe_draw(gg)
    # Ensure the canvas is rendered
    fig_tmp.set_dpi(dpi)
    fig_tmp.canvas.draw()

    # Copy the RGBA buffer into the target axes
    buf = np.asarray(fig_tmp.canvas.buffer_rgba())
    ax.imshow(buf)
    ax.axis("off")

    if close:
        plt.close(fig_tmp)


def gg_subplots(plots, ncols=2, figsize=(12, 6), dpi=200, titles=None, tight=True, gridspec_kw=None):
    """
    Arrange a list of plotnine ggplots in a Matplotlib grid.

    Parameters
    ----------
    plots   : list[ggplot]     plots to arrange
    ncols   : int              columns in the grid
    figsize : (w, h)           outer figure size in inches
    dpi     : int              DPI for crisp rendering of each ggplot
    titles  : list[str]|None   optional per-subplot titles
    tight   : bool             use tight_layout() at the end
    gridspec_kw : dict         additional kwargs for plt.subplots() gridspec

    Returns
    -------
    fig : matplotlib.figure.Figure
    axes : 2D ndarray of Axes (flattened if nrows==1 or ncols==1)
    """
    k = len(plots)
    nrows = math.ceil(k / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, dpi=dpi, gridspec_kw=gridspec_kw)
    axes = np.atleast_1d(axes).ravel()

    # Render each ggplot into its slot
    for i, gg in enumerate(plots):
        render_ggplot_to_ax(gg, axes[i], dpi=dpi)
        if titles and i < len(titles) and titles[i] is not None:
            axes[i].set_title(titles[i], pad=6)

    # Hide any unused slots
    for j in range(k, len(axes)):
        axes[j].set_visible(False)

    if tight:
        plt.tight_layout()
    return fig, axes


"""
LIMITATIONS:

- Subplots are rasters (images), not vector graphics. When saving as PDF/SVG, those panels won’t be vector elements.

- If you need fully vector output, either:
    - Use faceting inside a single ggplot (when feasible), or
    - Combine separate vector PDFs with a layout tool (LaTeX, Inkscape, Illustrator)

NOTES: 

- Legends are embedded inside each plot image. 
- A shared legend across subplots is non-trivial because we’re pasting images.

"""