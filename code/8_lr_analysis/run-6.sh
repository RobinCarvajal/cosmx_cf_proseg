#!/usr/bin/env bash

# Continue even if one niche/comparison fails
set +e

MAIN_DIR="/mnt/data/project0062/proseg_data"

cd "${MAIN_DIR}" || exit 1

SCRIPT_PATH="code/8_lr_analysis/6-results-comparison-final-cli.R"
CC_FUNCTIONS_PATH="${MAIN_DIR}/code/8_lr_analysis/cellchat-functions.R"

NICHES_DIR="${MAIN_DIR}/results/comb-full/lr-manual-niches"

# Get all niche directory names
mapfile -t NICHE_IDS < <(
  find "${NICHES_DIR}" \
    -mindepth 1 \
    -maxdepth 1 \
    -type d \
    -printf '%f\n'
)

# Define the three pairwise comparisons
# Format: REF CONT
COMPARISONS=(
  "CTRL CF"
  "CTRL CF_ETI"
  "CF CF_ETI"
)

for NICHE_ID in "${NICHE_IDS[@]}"; do

  echo "=================================================="
  echo "Processing niche: ${NICHE_ID}"
  echo "=================================================="

  INPUT_DIR="results/comb-full/lr-manual-niches/${NICHE_ID}"
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

      echo "Completed: ${NICHE_ID} — ${REF} vs ${CONT}"
      rm -f "${TMP_LOG}"

    else

      echo "Failed: ${NICHE_ID} — ${REF} vs ${CONT}"

      {
        echo "=================================================="
        echo "Failed niche: ${NICHE_ID}"
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
done