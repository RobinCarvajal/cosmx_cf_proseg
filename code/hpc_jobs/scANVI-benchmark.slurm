#!/bin/bash
#SBATCH --account=project0062
#SBATCH --job-name=scanvi-bench
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=32
#SBATCH --mem=80G
#SBATCH --time=08:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --mail-user=robin.carvajal@glasgow.ac.uk
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# ---- Conda env ----
module load apps/miniforge
eval "$(conda shell.bash hook)"
conda activate rapids_singlecell   # <-- change if your env is different

# (Optional) keep CPU threads sane for dataloaders / BLAS
export OMP_NUM_THREADS=2
export MKL_NUM_THREADS=2
export UMAP_DISABLE_JAX=1   # silence JAX GPU warning if you’re not using it

# ---- Paths ----
SCRIPT_PATH="/mnt/data/project0062/cosmx_gray/code/3_cell_typing/comb/4_benchmark.py"  # <- check filename (typo?)
WORKDIR="/mnt/data/project0062/cosmx_gray"
mkdir -p logs
cd "$WORKDIR"

# ---- Quick logging ----
echo "Host: $(hostname)"
echo "CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES"
nvidia-smi -L || true
which python

# ---- Run (use srun inside sbatch) ----
srun python "$SCRIPT_PATH"
