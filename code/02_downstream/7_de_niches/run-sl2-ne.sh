#!/usr/bin/env bash

# SOFTWARE PATHS
searchlight_dir="/mnt/data/project0062/cosmx_gray/software_ext/SearchLight2"
searchlight_py="${searchlight_dir}/software/Searchlight2.py" # Searchligth script
r_path="/mnt/data/project0062/.conda/envs/sl2/bin/Rscript" # R script

# DATA PATHS
bg_path="${searchlight_dir}/backgrounds/human_GRCh38.p13.cosmx.tsv" # background path
celltypes_dir="/mnt/data/project0062/cosmx_gray/results/comb/niches_de" # celltypes dir
ss_path="/mnt/data/project0062//cosmx_gray/data/sample_sheet/SL2_sample_sheet.tsv"

# GROUPS
cf="CF"
cf_eti="CF_ETI"
ctrl="CTRL"

# celltype list
# capture all immediate subdirectory names as an array
mapfile -t celltypes_list < <(find "$celltypes_dir" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort)

# print (optional)
printf "Found %d cell types:\n" "${#celltypes_list[@]}"
printf "  %s\n" "${celltypes_list[@]}"

# Print the array (optional)
echo "${celltypes_list[@]}"

for celltype in "${celltypes_list[@]}"; do
    echo "${celltype}"

    out_dir="${celltypes_dir}/${celltype}/SL2_results_ne"
    rm -rf "${out_dir}"
    mkdir -p "${out_dir}"

    # DE TABLES
    de_cf_ctrl="${celltypes_dir}/${celltype}/de_${cf}_vs_${ctrl}_${celltype}.tsv"
    de_cf_eti_ctrl="${celltypes_dir}/${celltype}/de_${cf_eti}_vs_${ctrl}_${celltype}.tsv"
    de_cf_eti_cf="${celltypes_dir}/${celltype}/de_${cf_eti}_vs_${cf}_${celltype}.tsv"

    if ! python3 "${searchlight_py}" --r path="${r_path}" --out path="${out_dir}" \
        --bg  file="${bg_path}" \
        --em  file="${celltypes_dir}/${celltype}/ne_${celltype}.tsv" \
        --ss  file="${ss_path}" \
        --de  file="${de_cf_ctrl}",numerator="${cf}",denominator="${ctrl}" \
        --de  file="${de_cf_eti_ctrl}",numerator="${cf_eti}",denominator="${ctrl}" \
        --de  file="${de_cf_eti_cf}",numerator="${cf_eti}",denominator="${cf}" \
        --mde name="multiple",numerator="${cf}"*denominator="${ctrl}",numerator="${cf_eti}"*denominator="${ctrl}",numerator="${cf_eti}"*denominator="${cf}" \
        --ora file="${searchlight_dir}/gene_set_databases/GO_human_BP.csv",type=GO_bp \
        --ura file="${searchlight_dir}/upstream_regulator_databases/trrust.human.tsv",type=TRRUST
    then
        echo "[FAIL] ${celltype} -> skipping (continuing)"
        # optional: remove partial output folder
        rm -rf "${out_dir}"
        continue
    fi
done