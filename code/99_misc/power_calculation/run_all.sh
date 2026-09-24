#!/usr/bin/env bash
set -euo pipefail

# Inputs that are shared across all cell types
base_ct_de="/Volumes/robin_work/cosmx_gray/results/comb/power"
ss_table="/Volumes/robin_work/cosmx_gray/data/sample_sheet/sample_sheet.tsv"
design_string='~ run + condition'

# Path to your R script
rscript_path="DESeq2.R"

# Loop over each celltype folder inside ct_de
for celltype_dir in "${base_ct_de}"/*/; do
  [[ -d "${celltype_dir}" ]] || continue

  celltype="$(basename "${celltype_dir%/}")"

  # Expected input/output paths for this cell type
  em_table="${celltype_dir}/pb_counts.csv"
  de_out_path="${celltype_dir}/de_table.tsv"
  ne_out_path="${celltype_dir}/vst_table.tsv"

  # Skip if the expected pseudobulk file isn't there
  if [[ ! -f "${em_table}" ]]; then
    echo "[SKIP] ${celltype}: missing ${em_table}"
    continue
  fi

  echo "[RUN ] ${celltype}"
  Rscript "${rscript_path}" \
    --em "${em_table}" \
    --ss "${ss_table}" \
    --design "${design_string}" \
    --de_out "${de_out_path}" \
    --ne_out "${ne_out_path}"
done
