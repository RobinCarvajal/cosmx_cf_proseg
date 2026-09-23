#!/usr/bin/env bash
set -euo pipefail

# run from project base dir (cosmx_gray)
niches_var="novae_domains_11"
base_ct_de="results/comb/de_${niches_var}"
ss_table="data/sample_sheet/sample_sheet.tsv"

ref_list=("CTRL" "CTRL" "CF")
query_list=("CF" "CF_ETI" "CF_ETI")

batch="run"
group="condition"

rscript_path="/mnt/data/project0062/cosmx_gray/code/7_de_niches/2-deseq2.R"

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

    em_table="${celltype_dir}/pb_${celltype}.csv"
    de_out_path="${celltype_dir}/de_${comparison}_${celltype}.tsv"
    ne_out_path="${celltype_dir}/ne_${celltype}.tsv"
    bc_out_path="${celltype_dir}/bc_${celltype}.tsv"

    if [[ ! -f "${em_table}" ]]; then
      echo "[SKIP] ${celltype}: missing ${em_table}"
      continue
    fi

    echo "[RUN] ${celltype} (${comparison})"

    if ! Rscript "${rscript_path}" \
      --em "${em_table}" \
      --ss "${ss_table}" \
      --batch "${batch}" \
      --group "${group}" \
      --ref "${ref}" \
      --query "${query}" \
      --de_out "${de_out_path}" \
      --ne_out "${ne_out_path}" \
      --bc_out "${bc_out_path}"
    then
      echo "[FAIL] ${celltype} (${comparison}) -> skipping (continuing)"
      continue
    fi

  done
done