#!/usr/bin/env bash

set -euo pipefail

PROJECT_ROOT="/users/rjc18y/storage/cosmx_cf_proseg"
cd "${PROJECT_ROOT}"

ANALYSIS_NAME="lr-inter-niche"

INPUT_PATH="data/comb-full/h5ad/comb-full-annotated.h5ad"
OUTPUT_PATH="results/comb-full/${ANALYSIS_NAME}/seurat-obj.qs"


spatial ccc adata-to-seurat --input_obj_path "${INPUT_PATH}" --output_obj_path "${OUTPUT_PATH}"

