#!/usr/bin/env bash

set -euo pipefail

PROJECT_ROOT="/users/rjc18y/storage/cosmx_cf_proseg"
cd "${PROJECT_ROOT}"

ANALYSIS_NAME="lr-intra-niche"

INPUT_PATH="data/comb-full/h5ad/comb-full-annotated.h5ad"
OUTPUT_PATH="results/comb-full/${ANALYSIS_NAME}"


spatial ccc adata-to-seurat-subsets \
 --input_obj_path "${INPUT_PATH}" \
 --output_dir "${OUTPUT_PATH}" \
 --niches_col "cellcharter_niches"

