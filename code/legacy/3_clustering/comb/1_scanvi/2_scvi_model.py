
# %% libraries 
import os
import scanpy as sc
from pathlib import Path
import numpy as np
import rapids_singlecell as rsc
import scvi

# %% setting paths
MAIN_DIR_NAME = "cosmx_gray"
MAIN_DIR = next(p for p in Path.cwd().parents if (p / MAIN_DIR_NAME).exists()) / MAIN_DIR_NAME
print(MAIN_DIR)
os.chdir(MAIN_DIR)

# %% Seed value
SEED_VALUE = 42
# set NumPy RNG for consistency
np.random.seed(42)
scvi.settings.seed = SEED_VALUE

# %% object versions
CUR_OBJ_V = 'v1'
CUR_OBJ_PATH = MAIN_DIR / 'data' / 'comb' / 'h5ad' / f'comb_{CUR_OBJ_V}.h5ad'

# %% Load the unintegrated object
comb = sc.read_h5ad(CUR_OBJ_PATH)

# %% Training model
scvi.model.SCVI.setup_anndata(
    comb,
    batch_key="run", # wild variation between runs
    layer="counts",
    categorical_covariate_keys=[
        "sample_name",  # sample_name has some variation
        'donor',       # donor has some variation
    ],
)

model = scvi.model.SCVI(comb)
model

# %%
print(comb.obs['_scvi_labels'].value_counts())
print(comb.obs['_scvi_batch'].value_counts())

# %%
scvi.settings.dl_num_workers = 63
scvi.settings.dl_persistent_workers = True

# %%
model.train(
    max_epochs=100, # number of iterations (100 works fine)
    accelerator='gpu', # 'cpu' or 'mps'/'gpu'
    devices=4, # if using multiple gpus, could use -1 to use all
    strategy="ddp_find_unused_parameters_true", # non-interactive

)

# %% [markdown]
# Saving the model

# %%
# saving the model

from pathlib import Path

models_dir = Path('data/comb')
models_dir.mkdir(parents=True, exist_ok=True)

print(models_dir)

# %%
# saving the model
scvi_model_path = MAIN_DIR / 'data/comb' / 'scvi_model'
model.save(scvi_model_path, overwrite=True)


