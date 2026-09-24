#!/usr/bin/env bash

# Exit immediately if a command fails
set -e

# The purpose of this script is to run script 4
# for each niche object created previously.

# Activate conda environment
#source /data/miniforge3/etc/profile.d/conda.sh
#conda activate cellchat-env


MAIN_DIR="/mnt/data/project0062/proseg_data"
cd "${MAIN_DIR}"

ANALYSIS_TYPE="whole"

SCRIPT_PATH="code/8_lr_analysis/${ANALYSIS_TYPE}/4-create-cellchat-objs-cli.R"
CC_FUNCTIONS_PATH="code/8_lr_analysis/${ANALYSIS_TYPE}/cellchat-functions.R"


INPUT_OBJ_PATH="data/comb-full/rds/comb-full-annotated-cellchat.rds"
OUTPUT_DIR="results/comb-full/lr-${ANALYSIS_TYPE}"
mkdir -p "${OUTPUT_DIR}"

LABELS_COL="ann_lvl_3_refined"

Rscript "${SCRIPT_PATH}" \
  --main_dir "${MAIN_DIR}" \
  --cc_functions_path "${CC_FUNCTIONS_PATH}" \
  --input_obj_path "${INPUT_OBJ_PATH}" \
  --output_dir "${OUTPUT_DIR}" \
  --labels_col "${LABELS_COL}"

