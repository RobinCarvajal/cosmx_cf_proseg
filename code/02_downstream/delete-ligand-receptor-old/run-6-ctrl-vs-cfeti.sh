#!/usr/bin/env bash

# Continue even if one niche fails
set +e

# Activate conda environment
#source /data/miniforge3/etc/profile.d/conda.sh
#conda activate cellchat-env

# Groups
REF="CTRL"
CONT="CF_ETI"

MAIN_DIR="/mnt/data/project0062/cosmx_gray"

cd "${MAIN_DIR}" || exit 1

SCRIPT_PATH="code/8_lr_analysis/6-results-comparison-final-cli.R"
CC_FUNCTIONS_PATH="${MAIN_DIR}/code/8_lr_analysis/cellchat-functions.R"

NICHES_DIR="${MAIN_DIR}/results/comb/lr-manual-niches"

NICHE_IDS=$(find "${NICHES_DIR}" \
  -mindepth 1 \
  -maxdepth 1 \
  -type d \
  | xargs -n 1 basename)

echo "${NICHE_IDS}"

for NICHE_ID in ${NICHE_IDS}; do

  echo "Processing niche: ${NICHE_ID}"

  INPUT_DIR="results/comb/lr-manual-niches/${NICHE_ID}"
  INPUT_DIR_FULL="${MAIN_DIR}/${INPUT_DIR}"

  OUTPUT_DIR="results/comb/lr-manual-niches/${NICHE_ID}/${REF}_vs_${CONT}"

  FAILED_LOG="${INPUT_DIR_FULL}/failed.log"

  # Remove previous failed log if it exists
  rm -f "${FAILED_LOG}"

  # Temporary log capturing full R output
  TMP_LOG=$(mktemp)

  if Rscript "${SCRIPT_PATH}" \
    --main_dir "${MAIN_DIR}" \
    --cc_functions_path "${CC_FUNCTIONS_PATH}" \
    --input_dir "${INPUT_DIR}" \
    --output_dir "${OUTPUT_DIR}" \
    --ref "${REF}" \
    --cont "${CONT}" \
    > "${TMP_LOG}" 2>&1; then

    echo "Niche ${NICHE_ID} completed successfully"

    rm -f "${TMP_LOG}"

  else

    echo "Niche ${NICHE_ID} failed"

    {
      echo "=================================================="
      echo "Failed niche: ${NICHE_ID}"
      echo "Date: $(date)"
      echo "=================================================="
      echo ""
      cat "${TMP_LOG}"
      echo ""
    } >> "${FAILED_LOG}"

    rm -f "${TMP_LOG}"

  fi

done