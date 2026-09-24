
#### DOES NOT WORK YET ####

# %% Libraries
import os
import scanpy as sc
import scvi
import numpy as np
from pathlib import Path
import rapids_singlecell as rsc
from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection
import datetime as dt
import matplotlib.pyplot as plt

# %% Setting Paths
# setting paths
MAIN_DIR_NAME = "cosmx_gray"
MAIN_DIR = next(p for p in Path.cwd().parents if (p / MAIN_DIR_NAME).exists()) / MAIN_DIR_NAME
os.chdir(MAIN_DIR)

# %% Setting Seed
SEED_VALUE = 42
# set NumPy RNG for consistency
np.random.seed(42)

# %% load comb object

# --- config
SCVI_LATENT_KEY = "X_scVI"
SCANVI_LATENT_KEY = "X_scANVI"
timestamp = dt.datetime.now().strftime("%Y%m%d-%H%M%S")
out_dir = Path("results/benchmarks") / f"scanvi_benchmark_{timestamp}"
out_dir.mkdir(parents=True, exist_ok=True)

# --- load
comb = sc.read_h5ad("data/comb/h5ad/comb_v5.h5ad")

# --- run benchmark
bm = Benchmarker(
    comb,
    batch_key="slide_name",
    label_key="scANVI_ct",
    embedding_obsm_keys=["X_pca", SCVI_LATENT_KEY, SCANVI_LATENT_KEY],
    n_jobs=-1,
)
bm.benchmark()

# --- plot results table and save
# If your Benchmarker supports return_fig=True, use that; otherwise grab the current fig.
try:
    fig = bm.plot_results_table(min_max_scale=False, return_fig=True)
except TypeError:
    bm.plot_results_table(min_max_scale=False)
    fig = plt.gcf()

for ext in ("png", "svg", "pdf"):
    fig.savefig(out_dir / f"results_table.{ext}", dpi=300, bbox_inches="tight")
plt.close(fig)
