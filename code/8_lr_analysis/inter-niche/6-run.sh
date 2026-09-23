#!/usr/bin/env bash

# Continue even if one niche/comparison fails
set +e

MAIN_DIR="/mnt/data/project0062/proseg_data"

cd "${MAIN_DIR}" || exit 1

SCRIPT_PATH="code/8_lr_analysis/inter-niche/6-results-comparison-final-cli.R"
CC_FUNCTIONS_PATH="${MAIN_DIR}/code/8_lr_analysis/inter-niche/cellchat-functions.R"

# Define the three pairwise comparisons
# Format: REF CONT
COMPARISONS=(
  "CTRL CF"
  "CTRL CF_ETI"
  "CF CF_ETI"
)

INPUT_DIR="results/comb-full/lr-inter-niche"
INPUT_DIR_FULL="${MAIN_DIR}/${INPUT_DIR}"

for COMPARISON in "${COMPARISONS[@]}"; do

  read -r REF CONT <<< "${COMPARISON}"

  echo "Running comparison: ${REF} vs ${CONT}"

  OUTPUT_DIR="${INPUT_DIR}/${REF}_vs_${CONT}"
  OUTPUT_DIR_FULL="${MAIN_DIR}/${OUTPUT_DIR}"

  mkdir -p "${OUTPUT_DIR_FULL}"

  # Separate failure log for each comparison
  FAILED_LOG="${OUTPUT_DIR_FULL}/failed.log"
  rm -f "${FAILED_LOG}"

  # Temporary file capturing the complete R output
  TMP_LOG=$(mktemp)

  if Rscript "${SCRIPT_PATH}" \
    --main_dir "${MAIN_DIR}" \
    --cc_functions_path "${CC_FUNCTIONS_PATH}" \
    --input_dir "${INPUT_DIR}" \
    --output_dir "${OUTPUT_DIR}" \
    --ref "${REF}" \
    --cont "${CONT}" \
    > "${TMP_LOG}" 2>&1; then

    echo "Completed: ${REF} vs ${CONT}"
    rm -f "${TMP_LOG}"

  else

    echo "Failed: ${REF} vs ${CONT}"

    {
      echo "=================================================="
      echo "Reference: ${REF}"
      echo "Contrast: ${CONT}"
      echo "Date: $(date)"
      echo "=================================================="
      echo ""
      cat "${TMP_LOG}"
      echo ""
    } > "${FAILED_LOG}"

    rm -f "${TMP_LOG}"

  fi

done
