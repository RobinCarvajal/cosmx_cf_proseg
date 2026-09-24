#!/usr/bin/env bash
set -euo pipefail

# this needs to run from the project base dir (cosmx_gray)

# Inputs that are shared across all cell types
base_ct_de="results/comb/ct_de"
ss_table="data/sample_sheet/sample_sheet.tsv"
ref_list=("CTRL" "CTRL" "CF")
query_list=("CF" "CF_ETI" "CF_ETI")
batch="run"
group="condition"


comparison="${query}_vs_${ref}"

# Path to your R script
rscript_path="/mnt/data/project0062/cosmx_gray/code/5_de_celltype/2-deseq2.R"

# Loop over each celltype folder inside ct_de
for celltype_dir in "${base_ct_de}"/*/; do
  [[ -d "${celltype_dir}" ]] || continue

  celltype="$(basename "${celltype_dir%/}")"

  # Expected input/output paths for this cell type
  em_table="${celltype_dir}/pb_${celltype}.csv"
  de_out_path="${celltype_dir}/de_${comparison}_${celltype}.tsv"
  ne_out_path="${celltype_dir}/ne_${celltype}.tsv"
  bc_out_path="${celltype_dir}/bc_${celltype}.tsv"

  # Skip if the expected pseudobulk file isn't there
  if [[ ! -f "${em_table}" ]]; then
    echo "[SKIP] ${celltype}: missing ${em_table}"
    continue
  fi

  echo "[RUN] ${celltype}"
  Rscript "${rscript_path}" \
    --em "${em_table}" \
    --ss "${ss_table}" \
    --batch "${batch}" \
    --group "${group}" \
    --ref "${ref}" \
    --query "${query}" \
    --de_out "${de_out_path}" \
    --ne_out "${ne_out_path}" \
    --bc_out "${bc_out_path}"
done
