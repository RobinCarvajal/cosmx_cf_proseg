#!/usr/bin/env bash
set -euo pipefail

# run from project base dir (cosmx_gray)

cd "/mnt/data/project0062/proseg_data"

base_ct_de="results/comb-full/de_ct"
ss_table="data/sample_sheet/sample_sheet.tsv"

ref_list=("CTRL" "CTRL" "CF")
query_list=("CF" "CF_ETI" "CF_ETI")

batch="run"
group="condition"

rscript_path="code/5_de_celltypes/2-deseq2.R"

# sanity check lengths match
if [[ "${#ref_list[@]}" -ne "${#query_list[@]}" ]]; then
  echo "ERROR: ref_list and query_list must have the same length" >&2
  exit 1
fi

# loop over comparisons (by index)
for i in "${!ref_list[@]}"; do
  ref="${ref_list[$i]}"
  query="${query_list[$i]}"
  comparison="${query}_vs_${ref}"

  echo "=== Comparison: ${comparison} ==="

  # loop over each celltype folder
  for celltype_dir in "${base_ct_de}"/*/; do
    [[ -d "${celltype_dir}" ]] || continue
    celltype="$(basename "${celltype_dir%/}")"

    # if [[ "${celltype}" == "hillock" ]]; then # skip hillock
    #   echo "[SKIP] hillock"
    #   continue
    # fi

    em_table="${celltype_dir}/pb_${celltype}.csv"
    de_out_path="${celltype_dir}/de_${comparison}_${celltype}.tsv"
    ne_out_path="${celltype_dir}/ne_${celltype}.tsv"
    bc_out_path="${celltype_dir}/bc_${celltype}.tsv"

    if [[ ! -f "${em_table}" ]]; then
      echo "[SKIP] ${celltype}: missing ${em_table}"
      continue
    fi

    echo "[RUN] ${celltype} (${comparison})"
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
done