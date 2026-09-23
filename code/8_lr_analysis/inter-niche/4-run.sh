#!/usr/bin/env bash

# Exit immediately if a command fails
set -e

# The purpose of this script is to run script 4
# for each niche object created previously.

# Activate conda environment
#source /data/miniforge3/etc/profile.d/conda.sh
#conda activate cellchat-env


MAIN_DIR="/mnt/data/project0062/proseg_data"

SCRIPT_PATH="${MAIN_DIR}/code/8_lr_analysis/inter-niche/4-create-cellchat-objs-cli.R"
CC_FUNCTIONS_PATH="${MAIN_DIR}/code/8_lr_analysis/inter-niche/cellchat-functions.R"


INPUT_OBJ_PATH="data/comb-full/rds/comb-full-annotated-cellchat.rds"
OUTPUT_DIR="results/comb-full/lr-inter-niche"
mkdir -p "${OUTPUT_DIR}"

Rscript "${SCRIPT_PATH}" \
  --main_dir "${MAIN_DIR}" \
  --cc_functions_path "${CC_FUNCTIONS_PATH}" \
  --input_obj_path "${INPUT_OBJ_PATH}" \
  --output_dir "${OUTPUT_DIR}"

