# %% Libraries
import os
import time
from pathlib import Path

import scanpy as sc
import pertpy as pt

# Start timer
start_time = time.perf_counter()

# %% Check GPU availability
try:
    import cupy as cp
    GPU_AVAIL = cp.cuda.runtime.getDeviceCount() > 0
except (ImportError, RuntimeError):
    GPU_AVAIL = False

if GPU_AVAIL:
    import rapids_singlecell as rsc

print(f"GPU available: {GPU_AVAIL}")

# %% Setting Paths
MAIN_DIR_NAME = "cosmx_cf_proseg"
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

comb = sc.read_h5ad(COMB_PATH)

# %% Milo
milo = pt.tl.Milo()
mdata = milo.load(comb)

if GPU_AVAIL:
    rsc.pp.neighbors(
        mdata["rna"],
        use_rep="X_resolvi",
        n_neighbors=150,
    )
else:
    sc.pp.neighbors(
        mdata["rna"],
        use_rep="X_resolvi",
        n_neighbors=150,
    )

milo.make_nhoods(
    mdata["rna"],
    prop=0.1,
)

mdata["rna"].obsm["nhoods"]

# %% Save MuData
out_dir = MAIN_DIR / "data" / "comb" / "h5mu" / "comb.h5mu"
out_dir.parent.mkdir(parents=True, exist_ok=True)
mdata.write(out_dir)

# %% Runtime
elapsed = time.perf_counter() - start_time

hours = int(elapsed // 3600)
minutes = int((elapsed % 3600) // 60)
seconds = elapsed % 60

runtime = (
    f"Runtime: {hours:02d}:{minutes:02d}:{seconds:05.2f}\n"
    f"Total seconds: {elapsed:.2f}\n"
    f"GPU available: {GPU_AVAIL}\n"
)

print(runtime)

runtime_file = out_dir.parent / "comb_milo_runtime.txt"
runtime_file.write_text(runtime)