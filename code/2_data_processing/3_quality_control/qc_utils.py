
# %% library imports
import scanpy as sc
import squidpy as sq
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

import warnings
from anndata._core.anndata import ImplicitModificationWarning

warnings.filterwarnings("ignore", category=ImplicitModificationWarning)


# %% flag cell counts

def qc_cell_counts(
    obj, 
    cell_cts_thr: int = 20, 
    cell_cts_var: str = "nCount_RNA"
):
    """
    Flag cells based on a *minimum counts* threshold and record pass/fail.

    Parameters
    ----------
    obj : AnnData
        AnnData object whose `.obs` contains per-cell counts.
    cell_cts_thr : int, default=20
        Minimum number of counts required for a cell to *pass* QC.
    cell_cts_var : str, default="nCount_RNA"
        Column in `obj.obs` holding per-cell counts.

    Side effects
    ------------
    - Adds/overwrites `obj.obs["qc_cell_counts"]`:
        1) First set as boolean (`True` = passes threshold).
        2) Then converted to string labels: `"pass"` / `"fail"` (categorical).

    Returns
    -------
    AnnData
        The same object, modified in place (returned for convenience/chaining).
    """
    import numpy as np  # ensure np is available inside the function

    # Step 1: boolean flag — True if cell meets or exceeds the minimum threshold
    # (i.e., passes QC), False otherwise.
    obj.obs["qc_cell_counts"] = obj.obs[cell_cts_var] >= cell_cts_thr

    # Step 2: map boolean to human-readable labels and store as categorical.
    # True  -> "pass"
    # False -> "fail"
    obj.obs["qc_cell_counts"] = np.where(
        obj.obs["qc_cell_counts"], "pass", "fail"
    )

    return obj


# %% qc fov_counts

def qc_fov_counts(
    obj, 
    fov_cts_thr: float = 40, 
    fov_cts_var: str = "nCount",
    fov_ncells_var: str = "nCell",
):
    """
    QC by average counts-per-cell at the FOV level.

    Adds to `obj.obs`:
      - 'avg_fov_cts_per_cell': average counts per cell for the FOV
        (computed as `fov_cts_var / fov_ncells_var`)
      - 'qc_fov_counts_bool': boolean flag (True = PASS if avg >= threshold)
      - 'qc_fov_counts': 'pass'/'fail' categorical label

    Parameters
    ----------
    obj : AnnData
        AnnData object whose `.obs` contains FOV-level columns.
    fov_cts_thr : float, default=40
        Minimum average counts-per-cell (FOV-level) required to pass QC.
    fov_cts_var : str, default='nCount'
        Column in `.obs` holding the **FOV total counts** replicated per cell
        (or a precomputed per-FOV metric aligned to each row).
    fov_ncells_var : str, default='nCell'
        Column in `.obs` holding the **number of cells in the FOV** replicated
        per cell.

    Returns
    -------
    AnnData
        Same object, modified in place (returned for convenience/chaining).

    Notes
    -----
    This implementation assumes `fov_cts_var` and `fov_ncells_var` are FOV-level
    values repeated for every cell in that FOV. If instead you only have
    per-cell counts and want the *mean per FOV* attached to each row, compute
    it via a groupby/transform on an FOV ID column.
    """
    import numpy as np

    obs = obj.obs

    # Average counts per cell in the FOV (assumes FOV-level totals replicated per row)
    obs["avg_fov_cts_per_cell"] = obs[fov_cts_var] / obs[fov_ncells_var]

    # Boolean flag: PASS if the FOV's avg counts-per-cell meets/exceeds threshold
    obs["qc_fov_counts_bool"] = obs["avg_fov_cts_per_cell"] >= fov_cts_thr

    # Human-readable categorical label derived from the boolean
    obs["qc_fov_counts"] = np.where(
        obs["qc_fov_counts_bool"], "pass", "fail"
    )

    return obj


# %% complexity 

def qc_complexity(
        obj, 
        gene_var='nFeature_RNA', 
        cell_cts_var='nCount_RNA', 
        complexity_thr=1
):
    """
    qc cells based on library complexity.

    Library complexity is approximated as the ratio of total counts (trancripts) 
    to the number of detected genes (features). Cells with unusually high 
    counts per gene may reflect low-complexity libraries.

    Parameters
    ----------
    obj : AnnData
        AnnData object containing single-cell or spatial transcriptomics data.
    gene_var : str, optional (default='nFeature_RNA')
        Column in obj.obs giving the number of detected genes per cell.
    cell_cts_var : str, optional (default='nCount_RNA')
        Column in obj.obs giving the total number of counts/UMIs per cell.
    complexity_thr : float, optional (default=1)
        Maximum allowed counts-per-gene ratio for a cell to be flagged. 
        Cells with ratio < threshold are flagged as low complexity.

    Returns
    -------
    obj : AnnData
        The same AnnData object with two new columns in .obs:
        - 'cell_complexity': counts-per-gene ratio for each cell
        - 'qc_complexity': True if cell passes the threshold, False otherwise
    """

    # Compute library complexity per cell:
    # ratio = total counts / number of detected genes
    obj.obs['cell_complexity'] = obj.obs[cell_cts_var] / obj.obs[gene_var]

    # qc cells that have low complexity:
    # True if ratio <= threshold (kept), False if ratio > threshold (filtered out)
    obj.obs['qc_complexity'] = obj.obs['cell_complexity'] >= complexity_thr

    # set as categorical (pass/fail)
    obj.obs['qc_complexity'] = np.where(
        obj.obs['qc_complexity'], "pass", "fail"
    )

    return obj


# %% negative counts 

def qc_negprobes(
    obj, 
    neg_prefix="Neg",
    cell_cts_var="nCount_RNA", 
    negprobe_thr=0.1
):
    """
    Perform QC filtering of cells based on negative probe counts.

    Parameters
    ----------
    obj : AnnData
        AnnData object containing gene expression data.
    neg_prefix : str, default="Neg"
        Prefix used to identify negative control probes in `var_names`.
    cell_cts_var : str, default="nCount_RNA"
        Name of the column in `obj.obs` with per-cell total RNA counts.
    negprobe_thr : float, default=0.1
        Threshold fraction of negative probe counts relative to total RNA counts.
        Cells with fraction >= threshold are marked as "fail".

    Returns
    -------
    obj : AnnData
        Same object with additional QC metrics stored in `.obs`:
        - "negprobes_count": total counts from negative probes.
        - "negprobes_percent": fraction of negative probe counts relative to total RNA.
        - "qc_negprobes": categorical label ("pass"/"fail").
    """

    # Identify all negative probes (genes/features) whose names start with the prefix
    neg_probes = obj.var_names[obj.var_names.str.startswith(neg_prefix)]

    # Sum expression counts of negative probes for each cell
    neg_counts = obj[:, neg_probes].X.sum(axis=1)

    # Convert counts to a dense 1D numpy array (safe even if X is sparse)
    import numpy as np
    neg_counts = np.asarray(neg_counts).ravel()

    # Store per-cell negative probe counts in `.obs`
    obj.obs["negprobes_count"] = neg_counts

    # Calculate the fraction of negative probe counts relative to total RNA counts
    obj.obs["negprobes_percent"] = (
        obj.obs["negprobes_count"] / obj.obs[cell_cts_var]
    )

    # Replace any NaN values (e.g., division by zero for empty cells) with 0
    obj.obs.fillna({"negprobes_percent":0}, inplace=True)

    # Flag cells that pass/fail QC:
    # pass if fraction of negative probes < threshold, fail otherwise
    obj.obs["qc_negprobes"] = obj.obs["negprobes_percent"] < negprobe_thr

    # Convert boolean QC results into categorical labels ("pass" / "fail")
    obj.obs['qc_negprobes'] = np.where(
        obj.obs['qc_negprobes'], "pass", "fail"
    )

    return obj

# %%

def qc_borders(
    obj,
    fov_names_col: str = "fov_name",
    x_col: str = "CenterX_global_px",
    y_col: str = "CenterY_global_px",
    rel_border_width: float = 0.01
):
    """
    Flag cells that are located near the border of their Field of View (FOV).

    Parameters
    ----------
    obj : AnnData
        AnnData object containing per-cell spatial metadata.
    fov_names_col : str, default="fov_name"
        Column in `.obs` that specifies the FOV identity for each cell.
    x_col : str, default="CenterX_global_px"
        Column in `.obs` with the global X coordinate of the cell centroid.
    y_col : str, default="CenterY_global_px"
        Column in `.obs` with the global Y coordinate of the cell centroid.
    rel_border_width : float, default=0.01
        Relative width (fraction of FOV size) to define the border region.  
        Example: `0.01` = border is 1% of FOV width/height from each edge.

    Returns
    -------
    obj : AnnData
        The same AnnData object with a new column in `.obs`:
        - "qc_borders": categorical label ("pass"/"fail") for border QC.
    """

    # Compute per-FOV minimum and maximum X and Y values.
    # `transform("min")` and `transform("max")` align results back to each row,
    # so every cell gets the min/max of its own FOV.
    x_min = obj.obs.groupby(fov_names_col)[x_col].transform("min")
    x_max = obj.obs.groupby(fov_names_col)[x_col].transform("max")
    y_min = obj.obs.groupby(fov_names_col)[y_col].transform("min")
    y_max = obj.obs.groupby(fov_names_col)[y_col].transform("max")

    # Define the absolute border widths in pixels for each FOV
    # by multiplying the FOV span by the relative width factor.
    x_w = (x_max - x_min) * rel_border_width
    y_w = (y_max - y_min) * rel_border_width

    # Identify cells that fall inside the border region:
    # - Close to left/right border (x near min/max)
    # - Close to top/bottom border (y near min/max)
    is_border = (
        (obj.obs[x_col] <= (x_min + x_w)) |
        (obj.obs[x_col] >= (x_max - x_w)) |
        (obj.obs[y_col] <= (y_min + y_w)) |
        (obj.obs[y_col] >= (y_max - y_w))
    )

    # Store boolean mask in `.obs`
    obj.obs["qc_borders"] = is_border

    # Convert boolean mask into categorical pass/fail labels
    # Cells at the border = "fail", others = "pass"
    obj.obs["qc_borders"] = np.where(
        obj.obs["qc_borders"], "fail", "pass"
    )

    return obj

# %% qc pass

def qc_pass(obj):
    """
    Combine multiple QC criteria into a final overall QC pass/fail decision.

    Parameters
    ----------
    obj : AnnData
        AnnData object containing QC results stored in `.obs`.
        Must already include the following QC columns:
        - "qc_cell_counts"
        - "qc_fov_counts"
        - "qc_complexity"
        - "qc_negprobes"
        - "qc_borders"

    Returns
    -------
    obj : AnnData
        The same AnnData object with an additional column:
        - "qc_pass": categorical label ("pass"/"fail") representing
          whether a cell passed all QC checks.
    """

    import pandas as pd

    # Boolean QC: a cell passes only if it passes all individual QC filters.
    # Each qc_* column is expected to be "pass" or "fail".
    obj.obs['qc_pass'] = (
        (obj.obs['qc_cell_counts'] == 'pass') &
        (obj.obs['qc_fov_counts'] == 'pass') &
        (obj.obs['qc_complexity'] == 'pass') &
        (obj.obs['qc_negprobes'] == 'pass') &
        (obj.obs['qc_borders'] == 'pass')
    )

    # Convert boolean (True/False) into categorical "pass"/"fail"
    obj.obs['qc_pass'] = np.where(obj.obs['qc_pass'], 'pass', 'fail')

    return obj



# %% final function 
def run_qc(
        obj,
        # shared parameters
        gene_var='nFeature_RNA',
        cell_cts_var="nCount_RNA", 
        fov_cts_var="nCount", 
        fov_ncells_var='nCell',
        # cell count parameters
        cell_cts_thr=20,
        # fov count parameters
        fov_names_col='fov_name', 
        fov_cts_thr=40, 
        # complexity parameters
        complexity_thr=1, 
        # negative probe parameters
        neg_prefix="Neg",
        negprobe_thr=0.1, 
        # border cell parameters
        x_col='CenterX_global_px', 
        y_col='CenterY_global_px', 
        rel_border_width=0.01, 
):
    qc_cell_counts(
        obj, 
        cell_cts_thr=cell_cts_thr, 
        cell_cts_var=cell_cts_var
    )

    qc_fov_counts(
        obj,
        fov_cts_thr=fov_cts_thr, 
        fov_cts_var=fov_cts_var, 
        fov_ncells_var=fov_ncells_var
    )

    qc_complexity(
        obj, 
        gene_var=gene_var, 
        cell_cts_var=cell_cts_var, 
        complexity_thr=complexity_thr
    )

    qc_negprobes(
        obj, 
        neg_prefix=neg_prefix, 
        cell_cts_var=cell_cts_var, 
        negprobe_thr=negprobe_thr
    )

    qc_borders(
        obj, 
        fov_names_col=fov_names_col, 
        x_col=x_col, 
        y_col=y_col, 
        rel_border_width=rel_border_width
    )

    qc_pass(obj)
    
    return obj


# %% qc filter

def qc_filter(obj):
    """
    Filter the AnnData object to retain only cells and features that 
    passed all quality control (QC) checks.

    Parameters
    ----------
    obj : AnnData
        AnnData object containing a 'qc_pass' column in `.obs` with
        categorical labels ("pass"/"fail") for each cell.

    Returns
    -------
    filtered_obj : AnnData
        A new AnnData object with:
        - Only cells that passed QC retained.
        - Control probes ("System*") and negative probes ("Neg*") removed from features.
    """

    # Ensure the 'qc_pass' column exists before filtering
    if 'qc_pass' not in obj.obs.columns:
        raise ValueError("The AnnData object must contain a 'qc_pass' column in .obs.")

    # Keep only cells marked as 'pass'
    # `.copy()` ensures this is a standalone AnnData object (not a view).
    filtered_obj = obj[obj.obs['qc_pass'] == 'pass'].copy()

    # Identify features (genes/probes) that are system controls or negative probes
    control_probes = obj.var_names[
        obj.var_names.str.startswith('System') |  # system control probes
        obj.var_names.str.startswith('Neg')      # negative probes
    ]

    # Remove those control/negative probes from the dataset
    filtered_obj = filtered_obj[:, ~filtered_obj.var_names.isin(control_probes)].copy()

    return filtered_obj