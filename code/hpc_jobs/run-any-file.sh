#!/bin/bash
#SBATCH --account=project0062
#SBATCH --job-name=SL2 # change name too 
#SBATCH --partition=gpu # gpu/gpuplus
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

# ---- Conda env ---- #
module load apps/miniforge
eval "$(conda shell.bash hook)"
conda activate sl2

# ---- Paths ---- #
SCRIPT_PATH="/mnt/data/project0062/proseg_data/code/5_de_celltypes/run-sl2-ne.sh"
WORKDIR="/mnt/data/project0062/proseg_data"

cd "$WORKDIR"

# ---- Run the script ---- #
srun bash "$SCRIPT_PATH"
