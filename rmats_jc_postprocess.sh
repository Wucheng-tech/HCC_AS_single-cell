#!/usr/bin/env bash
set -euo pipefail


# bash /public/home/jiangzh/projects/hcc_single_cell/sc_rmats/rmats_jc_postprocess.sh /public/home/wucheng/Analysis/AS/scRNA/scFAST/Sample/celltype_split/merged/T_N1/Tumor_Normal

# ===============================
# Usage check
# ===============================
# if [ $# -ne 2 ]; then
#     echo "Usage: $0 <rmats_res_directory> <script_path>"
#     exit 1
# fi

rmats_res="$1"
# script_path="$2"

# ===============================
# Locate rMATS root directory
# ===============================
rmats_root_dir=$(dirname "$(readlink -f "$(which rmats.py)")")

# ===============================
# AS types to process
# ===============================
AS_TYPES=(SE RI MXE A3SS A5SS)

# ===============================
# Main loop
# ===============================
for as_type in "${AS_TYPES[@]}"; do
    echo ">>> Processing AS type: ${as_type}"

    # ostat_fdr="${rmats_res}/tmp/JC_${as_type}/rMATS_result_FDR.txt"
    # ostat_il="${rmats_res}/tmp/JC_${as_type}/rMATS_result_I-L.txt"
    new_file="${rmats_res}/tmp/JC_${as_type}/rMATS_result_new.txt"
    gtf_file="${rmats_res}/fromGTF.${as_type}.txt"
    out_file="${rmats_res}/${as_type}.MATS.JC.new.txt"

    # # -------------------------------
    # # File existence check
    # # -------------------------------
    # if [[ ! -s "$ostat_fdr" || ! -s "$ostat_il" || ! -s "$gtf_file" ]]; then
    #     echo "  [WARN] Missing input files for ${as_type}, skip."
    #     continue
    # fi

    # # -------------------------------
    # # 1. sort FDR by numeric ID + paste with I-L
    # # -------------------------------
    # paste \
    # <( {
    #     # header
    #     head -n 1 "$ostat_fdr"

    #     # numeric ID
    #     tail -n +2 "$ostat_fdr" \
    #     | awk -F'\t' '$1 ~ /^[0-9]+$/ {print}' \
    #     | sort -t $'\t' -k1,1n

    #     # non-numeric ID (keep original order)
    #     tail -n +2 "$ostat_fdr" \
    #     | awk -F'\t' '$1 !~ /^[0-9]+$/ {print}'
    # } ) \
    # "$ostat_il" > "$new_file"

    # # -------------------------------
    # # 2. remove lines with column != 12
    # # -------------------------------
    # awk -F'\t' '
    # NR==1 { print; next }
    # NF == 12 { print }
    # ' "$new_file" > "${new_file}.tmp" && mv "${new_file}.tmp" "$new_file"

    # -------------------------------
    # 1. process tmp file0
    # -------------------------------
    python /public/home/jiangzh/projects/hcc_single_cell/sc_rmats/rmats_jc_postprocess.py ${rmats_res} ${as_type}
    # -------------------------------
    # 2. join with fromGTF
    # -------------------------------
    python "${rmats_root_dir}/rMATS_P/joinFiles.py" \
        "$gtf_file" \
        "$new_file" \
        0 0 \
        "$out_file"

    echo " Done: ${out_file}"
done

echo "All AS types processed."
