#!/usr/bin/env bash

set -euo pipefail

module load apps/miniforge
conda activate cellchat-env

PROJECT_ROOT="/users/rjc18y/storage/cosmx_cf_proseg"
cd "${PROJECT_ROOT}"

ANALYSIS_NAME="lr-intra-niche"


# Get all niche directories automatically
# NICHES_DIR="${PROJECT_ROOT}/results/comb-full/${ANALYSIS_NAME}"
# NICHE_IDS=($(find "${NICHES_DIR}" \
#   -mindepth 1 \
#   -maxdepth 1 \
#   -type d \
#   -exec basename {} \;))

# printf '%s\n' "${NICHE_IDS[@]}"

# if any groups failed run with only those niches
NICHE_IDS=(
  "Luminal_Airway"
  "Neutrophil_Rich"
  "Plasma_Rich"
)
printf '%s\n' "${NICHE_IDS[@]}"



# Loop through each niche object and run the script
for NICHE_ID in "${NICHE_IDS[@]}"; do

  echo "Processing niche: ${NICHE_ID}"

  INPUT_OBJ_PATH="results/comb-full/${ANALYSIS_NAME}/${NICHE_ID}/niche-obj.qs"
  OUTPUT_DIR="results/comb-full/${ANALYSIS_NAME}/${NICHE_ID}"


  spatial ccc create-cellchat-objs \
    --input_obj_path "${INPUT_OBJ_PATH}" \
    --output_dir "${OUTPUT_DIR}" \
    --labels_col "ann_lvl_3_refined" \
    --condition_col "condition" \
    --sample_col "sample_name" \
    --future_plan "sequential" \
    --future_workers 4 \
    --future_max_gb 380

done
