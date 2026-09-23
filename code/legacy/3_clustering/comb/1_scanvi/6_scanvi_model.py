
# %% Libraries
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

# %% object versions
CUR_OBJ_V = 'v4'
CUR_OBJ_PATH = MAIN_DIR / 'data' / 'comb' / 'h5ad' / f'comb_{CUR_OBJ_V}.h5ad'

# %% Setting scVI parameters
scvi.settings.dl_num_workers = 63
scvi.settings.dl_persistent_workers = True
scvi.settings.seed = SEED_VALUE

# %% Load the unintegrated object
comb = sc.read_h5ad(CUR_OBJ_PATH)

# %% Training model

# load scvi model
scvi_model_path = MAIN_DIR / 'data/comb' / 'scvi_model'
model = scvi.model.SCVI.load(
    scvi_model_path, 
    adata=comb,
    accelerator='gpu',
)

# %% create scanvi model from scvi model
SCVI_CT_KEY = "ct_scVI"
scanvi_model = scvi.model.SCANVI.from_scvi_model(
    model,
    adata=comb,
    labels_key=SCVI_CT_KEY,
    unlabeled_category="unknown",
)

# %% train scanvi model

## only use devices = -1 and strategy = "ddp_find_unused_parameters_true" in interactive environments
## if the script fails try using less gpus or cpu only
## also modify the slurm script accordingly
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2   # one process per GPU
#SBATCH --gres=gpu:2  

scanvi_model.train(
    max_epochs=100, # maybe try 100
    accelerator='gpu',
    devices=1, # if using multiple gpus, could use -1 to use all
    #strategy="ddp_find_unused_parameters_true", # non-interactive
) 

# %% saving the model
models_dir = Path('data/comb')
models_dir.mkdir(parents=True, exist_ok=True)

print(models_dir)

scanvi_model_path = MAIN_DIR / 'data/comb' / 'scanvi_model'
scanvi_model.save(scanvi_model_path, overwrite=True)

# %%
