#!/bin/bash
#SBATCH --account=project0062
#SBATCH --job-name=LR
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

# ---- Conda env ---- #
module load apps/miniforge
eval "$(conda shell.bash hook)"
conda activate cellchat-env

# ---- Paths ---- #
WORKDIR="/mnt/data/project0062/cosmx_gray"
cd "$WORKDIR"


LR_SCRIPTS_DIR="/mnt/data/project0062/cosmx_gray/code/8_lr_analysis"

SCRIPT1_PATH="${LR_SCRIPTS_DIR}/1-cc-adata.py"
SCRIPT2_PATH="${LR_SCRIPTS_DIR}/2-cc-adata-to-rds.R"
SCRIPT3_PATH="${LR_SCRIPTS_DIR}/3-create-niche-objs.R"
SCRIPT4_PATH="${LR_SCRIPTS_DIR}/run-4.sh"
SCRIPT6_1_PATH="${LR_SCRIPTS_DIR}/run-6-cf-vs-cfeti.sh"
SCRIPT6_2_PATH="${LR_SCRIPTS_DIR}/run-6-ctrl-vs-cf.sh"
SCRIPT6_3_PATH="${LR_SCRIPTS_DIR}/run-6-ctrl-vs-cfeti.sh"

# ---- Run the script ---- #

conda activate rsc_25.12
srun python "$SCRIPT1_PATH"

conda activate cellchat-env
srun Rscript "$SCRIPT2_PATH"
srun Rscript "$SCRIPT3_PATH"
srun bash "$SCRIPT4_PATH"
srun bash "$SCRIPT6_1_PATH"
srun bash "$SCRIPT6_2_PATH"
srun bash "$SCRIPT6_3_PATH"
