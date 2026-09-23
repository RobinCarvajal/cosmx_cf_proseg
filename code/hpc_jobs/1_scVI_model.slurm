#!/bin/bash
#SBATCH --account=project0062
#SBATCH --job-name=scVI-model
#SBATCH --partition=gpuplus
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4   # one process per GPU
#SBATCH --gres=gpu:4            # request 4 GPUs
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=08:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --mail-user=robin.carvajal@glasgow.ac.uk
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# ---- Conda env ----
module load apps/miniforge
eval "$(conda shell.bash hook)"
conda activate rapids_singlecell

# ---- Environment sanity ----
export OMP_NUM_THREADS=2
export MKL_NUM_THREADS=2
export UMAP_DISABLE_JAX=1
export TORCH_NCCL_BLOCKING_WAIT=1     # optional, helps debugging hangs
export NCCL_DEBUG=INFO                # verbose NCCL output (optional)

# ---- Paths ----
SCRIPT_PATH="/mnt/data/project0062/cosmx_gray/code/2_clustering/comb/scripts/2_scvi_model.py"
WORKDIR="/mnt/data/project0062/cosmx_gray"
mkdir -p logs
cd "$WORKDIR"

# ---- Quick logging ----
echo "Host: $(hostname)"
echo "CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES"
nvidia-smi -L || true
which python

# ---- Run (multi-GPU, one process per GPU) ----
# Lightning automatically handles DDP if you specify devices=4 in Trainer()
# Using srun ensures proper rank assignment across GPUs
srun python "$SCRIPT_PATH"
