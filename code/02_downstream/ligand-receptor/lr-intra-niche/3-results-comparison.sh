#!/usr/bin/env bash

# Continue if one niche/comparison fails, but fail fast on real setup errors
set -u

PROJECT_ROOT="/mnt/data/project0062/cosmx_cf_proseg"
cd "${PROJECT_ROOT}" || exit 1

module load apps/miniforge
conda activate cellchat-env

# ── Guard 1: 'spatial' must be resolvable ─────────────────────────────
command -v spatial >/dev/null 2>&1 \
  || { echo "ERROR: 'spatial' not found — conda env 'cellchat-env' did not activate." >&2; exit 1; }

ANALYSIS_TYPE="lr-intra-niche"
NICHES_DIR="${PROJECT_ROOT}/results/comb-full/${ANALYSIS_TYPE}"

# ── Guard 2: niche dir must exist and contain at least one niche ─────
if [[ ! -d "${NICHES_DIR}" ]] \
   || [[ -z "$(find "${NICHES_DIR}" -mindepth 1 -maxdepth 1 -type d -print -quit)" ]]; then
  echo "ERROR: no niche directories found under ${NICHES_DIR}" >&2
  exit 1
fi

readarray -t NICHE_IDS < <(find "${NICHES_DIR}" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort)

printf '%s\n' "${NICHE_IDS[@]}"

COMPARISONS=(
  "CTRL CF"
  "CTRL CF_ETI"
  "CF CF_ETI"
)

LABELS_CONFIG="code/02_downstream/ligand-receptor/${ANALYSIS_TYPE}/labels_config.qs"

# ── Guard 3: labels_config must exist ─────────────────────────────────
[[ -f "${LABELS_CONFIG}" ]] \
  || { echo "ERROR: labels_config not found at ${LABELS_CONFIG}" >&2; exit 1; }

# ── Clean up temp files on early exit ─────────────────────────────────
TMP_LOG=""
cleanup() { [[ -n "${TMP_LOG}" ]] && [[ -f "${TMP_LOG}" ]] && rm -f "${TMP_LOG}"; }
trap cleanup EXIT

for NICHE_ID in "${NICHE_IDS[@]}"; do

  echo "=================================================="
  echo "Processing niche: ${NICHE_ID}"
  echo "=================================================="

  INPUT_DIR="results/comb-full/${ANALYSIS_TYPE}/${NICHE_ID}"

  for COMPARISON in "${COMPARISONS[@]}"; do

    read -r REF CONT <<< "${COMPARISON}"
    echo "Running comparison: ${REF} vs ${CONT}"

    OUTPUT_DIR="${INPUT_DIR}/${REF}_vs_${CONT}"
    mkdir -p "${OUTPUT_DIR}"

    FAILED_LOG="${OUTPUT_DIR}/failed.log"
    rm -f "${FAILED_LOG}"

    TMP_LOG=$(mktemp)

    if spatial ccc results-comparison \
        --input_dir "${INPUT_DIR}" \
        --output_dir "${OUTPUT_DIR}" \
        --ref "${REF}" \
        --cont "${CONT}" \
        --labels_config "${LABELS_CONFIG}" \
        --future_plan "sequential" \
        --future_workers 1 \
        --future_max_gb 380 \
      > "${TMP_LOG}" 2>&1; then

      # Optionally keep a success log:
      # mv "${TMP_LOG}" "${OUTPUT_DIR}/run.log"
      echo "Completed: ${NICHE_ID} — ${REF} vs ${CONT}"

    else

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
      echo "Failed: ${NICHE_ID} — ${REF} vs ${CONT}"

    fi

    TMP_LOG=""

  done
done