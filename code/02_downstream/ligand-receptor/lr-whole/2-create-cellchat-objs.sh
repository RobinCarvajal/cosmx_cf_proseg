#!/usr/bin/env bash

set -euo pipefail

PROJECT_ROOT="/users/rjc18y/storage/cosmx_cf_proseg"
cd "${PROJECT_ROOT}"

ANALYSIS_NAME="lr-whole"

INPUT_PATH="results/comb-full/${ANALYSIS_NAME}/seurat-obj.qs"
OUTPUT_PATH="results/comb-full/${ANALYSIS_NAME}"


spatial ccc create-cellchat-objs \
    --input_obj_path "${INPUT_PATH}" \
    --output_dir "${OUTPUT_PATH}" \
    --labels_col "ann_lvl_3_refined" \
    --condition_col "condition" \
    --sample_col "sample_name" \
    --future_plan "sequential" \
    --future_workers 4 \
    --future_max_gb 380
