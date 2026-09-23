#!/bin/bash
#SBATCH --account=project0062
#SBATCH --job-name=ResolVI
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gres=gpu:1            # request 4 GPUs
#SBATCH --ntasks-per-node=1   # one process per GPU
#SBATCH --cpus-per-task=64
#SBATCH --mem=120G
#SBATCH --time=10:00:00
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
ulimit -n 65535  # increase file descriptor limit
export OMP_NUM_THREADS=2
export MKL_NUM_THREADS=2
export UMAP_DISABLE_JAX=1
export TORCH_NCCL_BLOCKING_WAIT=1     # optional, helps debugging hangs
export NCCL_DEBUG=INFO                # verbose NCCL output (optional)

# ---- Paths ----
SCRIPT_PATH="/mnt/data/project0062/cosmx_gray/code/3_cell_typing/comb/4_resolvi_model.py"  # <- check filename (typo?)
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
