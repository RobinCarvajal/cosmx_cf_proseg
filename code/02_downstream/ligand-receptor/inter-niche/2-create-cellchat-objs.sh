#!/usr/bin/env bash

set -euo pipefail

PROJECT_ROOT="/users/rjc18y/storage/cosmx_cf_proseg"
cd "${PROJECT_ROOT}"

ANALYSIS_NAME="lr-inter-niche"

INPUT_PATH="results/comb-full/${ANALYSIS_NAME}/seurat-obj.qs"
OUTPUT_PATH="results/comb-full/${ANALYSIS_NAME}"


spatial ccc create-cellchat-objs \ 
    --input_obj_path "${INPUT_PATH}" \
    --output_dir "${OUTPUT_PATH}" \
    --labels_col "cellcharter_niches" \
    --condition_col "condition" \
    --sample_col "sample_name" \
    --future_plan "multisession" \
    --future_workers 4 \
    --future_max_gb 300
