#!/usr/bin/env bash

set -euo pipefail

module load apps/miniforge
conda activate cellchat-env

PROJECT_ROOT="/users/rjc18y/storage/cosmx_cf_proseg"
cd "${PROJECT_ROOT}"

ANALYSIS_NAME="lr-inter-niche"

INPUT_PATH="results/comb-full/${ANALYSIS_NAME}"
OUTPUT_PATH="results/comb-full/${ANALYSIS_NAME}"
LABELS_CONFIG="code/02_downstream/ligand-receptor/inter-niche/labels_config.qs"

# Format: REF CONT
COMPARISONS=(
    "CTRL CF"
    "CTRL CF_ETI"
    "CF CF_ETI"
)

for COMPARISON in "${COMPARISONS[@]}"; do
    read -r REF CONT <<< "${COMPARISON}"
    echo "Running comparison: ${REF} vs ${CONT}"

    OUTPUT_DIR="${OUTPUT_PATH}/${REF}_vs_${CONT}"
    mkdir -p "${OUTPUT_DIR}"

    FAILED_LOG="${OUTPUT_DIR}/failed.log"
    RUN_LOG="${OUTPUT_DIR}/run.log"
    rm -f "${FAILED_LOG}"

    if spatial ccc results-comparison \
        --input_dir "${INPUT_PATH}" \
        --output_dir "${OUTPUT_DIR}" \
        --ref "${REF}" \
        --cont "${CONT}" \
        --labels_config "${LABELS_CONFIG}" \
        --future_plan "sequential" \
        --future_workers 1 \
        --future_max_gb 380 \
        > "${RUN_LOG}" 2>&1; then

        echo "Completed: ${REF} vs ${CONT}"
    else
        EXIT_CODE=$?
        echo "Failed: ${REF} vs ${CONT} (exit ${EXIT_CODE})"

        {
            echo "Failed comparison: ${REF} vs ${CONT}"
            echo "Exit code: ${EXIT_CODE}"
            echo "Date: $(date)"
            echo ""
            cat "${RUN_LOG}"
        } > "${FAILED_LOG}"
    fi
done