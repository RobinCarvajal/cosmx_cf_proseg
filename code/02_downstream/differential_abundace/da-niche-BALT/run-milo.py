# %% Libraries
import os
import time
import scanpy as sc
from pathlib import Path
#import rapids_singlecell as rsc
import pertpy as pt
import matplotlib.pyplot as plt
import numpy as np

# SEED
SEED_VALUE=42

# Start timer
start_time = time.perf_counter()

# %% Setting Paths
MAIN_DIR_NAME = "proseg_data"
MAIN_DIR = next(
    p for p in Path.cwd().parents
    if (p / MAIN_DIR_NAME).exists()
) / MAIN_DIR_NAME

os.chdir(MAIN_DIR)

# %% Object versions
COMB_V = "cellcharter"
COMB_PATH = (
    MAIN_DIR / "data" / "comb" / "h5ad" /
    f"comb-{COMB_V}.h5ad"
)

# Load comb obj
comb = sc.read_h5ad(COMB_PATH)

# subset to airway_epithelium niche 
NICHE_NAME = 'BALT'
niche_mask = comb.obs['cellcharter_niches'] == NICHE_NAME
niche_obj = comb[niche_mask].copy()

# %% Milo
milo = pt.tl.Milo()

# Initialize object for Milo analysis
mdata = milo.load(niche_obj)

# rsc.pp.neighbors(
#     mdata["rna"],
#     use_rep="X_resolvi",
#     n_neighbors=150,
# )

sc.pp.neighbors(
    mdata["rna"],
    use_rep="X_resolvi",
    n_neighbors=150,
    random_state=SEED_VALUE,
)

milo.make_nhoods(
    mdata["rna"],
    prop=0.1,
)

mdata["rna"].obsm["nhoods"]

# %% Save MuData
out_dir = MAIN_DIR / "data" / 'niches' / NICHE_NAME / 'h5mu'
out_dir.mkdir(parents=True, exist_ok=True)
mdata.write(
  out_dir / "mdata.h5mu"
)

# %% Runtime
elapsed = time.perf_counter() - start_time

hours = int(elapsed // 3600)
minutes = int((elapsed % 3600) // 60)
seconds = elapsed % 60

runtime = (
    f"Runtime: {hours:02d}:{minutes:02d}:{seconds:05.2f}\n"
    f"Total seconds: {elapsed:.2f}\n"
)

print(runtime)

runtime_file = out_dir / "milo_runtime.txt"

with open(runtime_file, "w") as f:
    f.write(runtime)

# the rest of the steps are less computationally heavy