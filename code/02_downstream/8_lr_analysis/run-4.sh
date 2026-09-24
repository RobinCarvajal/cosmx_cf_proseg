#!/usr/bin/env bash

# Exit immediately if a command fails
set -e

# The purpose of this script is to run script 4
# for each niche object created previously.

# Activate conda environment
#source /data/miniforge3/etc/profile.d/conda.sh
#conda activate cellchat-env


MAIN_DIR="/mnt/data/project0062/proseg_data"

SCRIPT_PATH="${MAIN_DIR}/code/8_lr_analysis/4-create-cellchat-objs-cli.R"
CC_FUNCTIONS_PATH="${MAIN_DIR}/code/8_lr_analysis/cellchat-functions.R"

# Get all niche directories automatically
# NICHES_DIR="${MAIN_DIR}/results/comb-full/lr-manual-niches"
# NICHE_IDS=($(find "${NICHES_DIR}" \
#   -mindepth 1 \
#   -maxdepth 1 \
#   -type d \
#   -exec basename {} \;))
#
# printf '%s\n' "${NICHE_IDS[@]}"

# For now, run only selected niches
NICHE_IDS=(
  "neutrophil_rich"
  "plasma_rich"
)

# Loop through each niche object and run the script
for NICHE_ID in "${NICHE_IDS[@]}"; do
  echo "Processing niche: ${NICHE_ID}"

  INPUT_OBJ_PATH="results/comb-full/lr-manual-niches/${NICHE_ID}/niche-obj.rds"
  OUTPUT_DIR="results/comb-full/lr-manual-niches/${NICHE_ID}"

  Rscript "${SCRIPT_PATH}" \
    --main_dir "${MAIN_DIR}" \
    --cc_functions_path "${CC_FUNCTIONS_PATH}" \
    --input_obj_path "${INPUT_OBJ_PATH}" \
    --output_dir "${OUTPUT_DIR}"
done
