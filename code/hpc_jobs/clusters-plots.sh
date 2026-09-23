#!/bin/bash
#SBATCH --account=project0062
#SBATCH --job-name=Clusters
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
conda activate rsc_25.12

# ---- Paths ----
SCRIPT_CLUSTERS="/mnt/data/project0062/cosmx_gray/code/3_clusters/comb/3_leiden.py"  # <- check filename (typo?)
SCRIPT_PLOTS="/mnt/data/project0062/cosmx_gray/results/comb/plots/clusters/code.py"
WORKDIR="/mnt/data/project0062/cosmx_gray"
mkdir -p logs

# ---- Run (multi-GPU, one process per GPU) ----

#srun python "$SCRIPT_CLUSTERS"
srun python "$SCRIPT_PLOTS"
