#!/bin/bash

## ==============================
## rmats analysis
## 1. Bulk HCC comparison: Tumor vs Normal
## ==============================

# Output files for BAM lists
tumor_out="/public/home/wucheng/Analysis/AS/bulk_seq/comp1/Tumor_BAM_list.txt"
normal_out="/public/home/wucheng/Analysis/AS/bulk_seq/comp1/Normal_BAM_list.txt"

# List of sample IDs
samples=("N-254609" "N-255952" "N-256677" "T-250891" "T-254165" "T-254609" "T-254959" "T-255350" "T-255639" "T-256677"
         "M-262103" "N-262103" "T-247958" "T-261911" "T-261948" "T-262103" "T-262253" "T-262303" "T-262343" "T-263241"
         "T-263307" "T-263375" "T-263844" "N-256287" "N-257590" "N-257789" "N-258192" "N-258490" "N-260699" "N-260865"
         "N-261176" "T-256287" "T-257590" "T-257789" "T-258192" "T-258490" "T-259195" "T-260699" "T-260865" "T-261176"
         "N-265910" "T-156963" "T-233586" "T-257590.2" "T-263401" "T-265288" "T-265910" "T-266731" "T-267289" "T-267569")

# Generate BAM file lists for tumor and normal
for i in "${samples[@]}"
do
    bam_path="/public/home/wucheng/Analysis/AS/bulk_seq/${i}/STAR/${i}_Aligned.sortedByCoord.out.bam"

    if [[ "$i" == T-* || "$i" == M-* ]]; then
        echo -n "$bam_path," >> "$tumor_out"
    elif [[ "$i" == N-* ]]; then
        echo -n "$bam_path," >> "$normal_out"
    fi
done

# Remove trailing comma
sed -i 's/,\s*$//' "$tumor_out"
sed -i 's/,\s*$//' "$normal_out"

# Run rmats for Bulk HCC Tumor vs Normal
compare_group="Tumor_Normal"
/public/home/jiangzh/miniconda3/envs/rmats/bin/rmats.py \
    --b1 "$tumor_out" \
    --b2 "$normal_out" \
    --gtf /public/home/wucheng/software/RNAediting/reference/gencode.v44.primary_assembly.annotation.gtf \
    --od /public/home/wucheng/Analysis/AS/bulk_seq/comp1/${compare_group} \
    --tmp /public/home/wucheng/Analysis/AS/bulk_seq/comp1/${compare_group}_tmp_rmats \
    -t paired \
    --variable-read-length \
    --readLength 150 \
    --nthread 20

####
rmats_filter=/public/home/wucheng/Analysis/AS/bulk_seq/software/rmats-turbo-tutorial-main/scripts/rmats_filtering.py
declare -a event_array=("SE" "MXE" "RI" "A5SS" "A3SS")
cd /public/home/wucheng/Analysis/AS/bulk_seq/comp1/Fit
for event in "${event_array[@]}"; do
    python $rmats_filter /public/home/wucheng/Analysis/AS/bulk_seq/comp/Tumor_Normal/${event}.MATS.JC.txt 10,0.05,0.95,0.01,0.05,0.5,0.2
done


#!/bin/bash

## =======================================================
## Bulk HCC comparison: Tumor vs Normal
## rmats analysis, filtering, and PSI matrix construction
## =======================================================

# -------------------------------
# 1. Define output BAM list files
# -------------------------------
tumor_out="/public/home/wucheng/Analysis/AS/bulk_seq/comp1/Tumor_BAM_list.txt"
normal_out="/public/home/wucheng/Analysis/AS/bulk_seq/comp1/Normal_BAM_list.txt"

# -------------------------------
# 2. Sample IDs
# -------------------------------
samples=("N-254609" "N-255952" "N-256677" "T-250891" "T-254165" "T-254609" "T-254959" "T-255350" "T-255639" "T-256677"
         "M-262103" "N-262103" "T-247958" "T-261911" "T-261948" "T-262103" "T-262253" "T-262303" "T-262343" "T-263241"
         "T-263307" "T-263375" "T-263844" "N-256287" "N-257590" "N-257789" "N-258192" "N-258490" "N-260699" "N-260865"
         "N-261176" "T-256287" "T-257590" "T-257789" "T-258192" "T-258490" "T-259195" "T-260699" "T-260865" "T-261176"
         "N-265910" "T-156963" "T-233586" "T-257590.2" "T-263401" "T-265288" "T-265910" "T-266731" "T-267289" "T-267569")

# -------------------------------
# 3. Generate BAM file lists
# -------------------------------
for i in "${samples[@]}"; do
    bam_path="/public/home/wucheng/Analysis/AS/bulk_seq/${i}/STAR/${i}_Aligned.sortedByCoord.out.bam"

    if [[ "$i" == T-* || "$i" == M-* ]]; then
        echo -n "$bam_path," >> "$tumor_out"
    elif [[ "$i" == N-* ]]; then
        echo -n "$bam_path," >> "$normal_out"
    fi
done

# Remove trailing commas
sed -i 's/,\s*$//' "$tumor_out"
sed -i 's/,\s*$//' "$normal_out"

# -------------------------------
# 4. Run rmats
# -------------------------------
compare_group="Tumor_Normal"
/public/home/jiangzh/miniconda3/envs/rmats/bin/rmats.py \
    --b1 "$tumor_out" \
    --b2 "$normal_out" \
    --gtf /public/home/wucheng/software/RNAediting/reference/gencode.v44.primary_assembly.annotation.gtf \
    --od /public/home/wucheng/Analysis/AS/bulk_seq/comp1/${compare_group} \
    --tmp /public/home/wucheng/Analysis/AS/bulk_seq/comp1/${compare_group}_tmp_rmats \
    -t paired \
    --variable-read-length \
    --readLength 150 \
    --nthread 20

# -------------------------------
# 5. Filter rmats events
# -------------------------------
rmats_filter=/public/home/wucheng/Analysis/AS/bulk_seq/software/rmats-turbo-tutorial-main/scripts/rmats_filtering.py
declare -a event_array=("SE" "MXE" "RI" "A5SS" "A3SS")

cd /public/home/wucheng/Analysis/AS/bulk_seq/comp1/Fit

for event in "${event_array[@]}"; do
    python $rmats_filter /public/home/wucheng/Analysis/AS/bulk_seq/comp1/${compare_group}/${event}.MATS.JC.txt 10,0.05,0.95,0.01,0.05,0.5,0.2
done

# -------------------------------
# 6. Merge all filtered events into IncLevel summary
# -------------------------------
input_dir="/public/home/wucheng/Analysis/AS/bulk_seq/comp1/Fit"
output_file="${input_dir}/IncLevel_summary.txt"
> "$output_file"  # clear output

for event in SE RI MXE A3SS A5SS; do
    input_file="${input_dir}/filtered_${event}.MATS.JC.txt"
    if [[ ! -f "$input_file" ]]; then
        echo "❌ File not found: $input_file"
        continue
    fi
    echo "✅ Processing: $input_file"

    # Find IncLevel1 and IncLevel2 columns
    header=$(head -n 1 "$input_file")
    IFS=$'\t' read -r -a cols <<< "$header"
    inc1_idx=0
    inc2_idx=0
    for i in "${!cols[@]}"; do
        if [[ "${cols[$i]}" == "IncLevel1" ]]; then
            inc1_idx=$((i+1))
        elif [[ "${cols[$i]}" == "IncLevel2" ]]; then
            inc2_idx=$((i+1))
        fi
    done

    if (( inc1_idx == 0 || inc2_idx == 0 )); then
        echo "⚠️ Skipping $event: IncLevel columns not found"
        continue
    fi

    # Extract IncLevel columns and add event ID
    awk -v e="$event" -v i1="$inc1_idx" -v i2="$inc2_idx" 'BEGIN{FS=OFS="\t"; n=0}
    NR>1 { n++; print e"_"n, $i1, $i2 }' "$input_file" >> "$output_file"
done
echo "🎉 Combined IncLevel summary created: $output_file"

# -------------------------------
# 7. Generate MetaInfo summary
# -------------------------------
output_meta="${input_dir}/MetaInfo_summary.txt"
> "$output_meta"

for event in SE RI MXE A3SS A5SS; do
    input_file="${input_dir}/filtered_${event}.MATS.JC.txt"
    if [[ ! -f "$input_file" ]]; then
        echo "❌ File not found: $input_file"
        continue
    fi
    echo "✅ Processing: $input_file"

    awk -v e="$event" 'BEGIN{FS=OFS="\t"; n=0}
    NR==1 { next }
    { n++; print e"_"n, $1, $2, $3, $4, $5 }' "$input_file" >> "$output_meta"
done
echo "🎉 MetaInfo summary created: $output_meta"

# -------------------------------
# 8. Generate PSI matrix
# -------------------------------
input_file="${input_dir}/IncLevel_summary.txt"
output_file="${input_dir}/IncLevel_summary_annotated.txt"
tumor_file="$tumor_out"
normal_file="$normal_out"

# Parse tumor sample IDs
tumor_samples=()
IFS=',' read -ra tumor_paths <<< "$(cat "$tumor_file")"
for path in "${tumor_paths[@]}"; do
    rel_path="${path##*/}"  # just file name
    id="${rel_path%.*}"     # remove extension
    tumor_samples+=("$id")
done

# Parse normal sample IDs
normal_samples=()
IFS=',' read -ra normal_paths <<< "$(cat "$normal_file")"
for path in "${normal_paths[@]}"; do
    rel_path="${path##*/}"
    id="${rel_path%.*}"
    normal_samples+=("$id")
done

# Convert to comma-separated strings
ntumor=${#tumor_samples[@]}
nnormal=${#normal_samples[@]}
tumor_str=$(IFS=,; echo "${tumor_samples[*]}")
normal_str=$(IFS=,; echo "${normal_samples[*]}")

# Annotate PSI matrix
awk -v tumors="$tumor_str" -v normals="$normal_str" -v ntumor="$ntumor" -v nnormal="$nnormal" '
BEGIN {
    FS=OFS="\t";
    split(tumors, tumor_arr, ",");
    split(normals, normal_arr, ",");
    printf "ID";
    for (i=1;i<=ntumor;i++) printf "\tPSI_Tumor_%s", tumor_arr[i];
    for (i=1;i<=nnormal;i++) printf "\tPSI_Normal_%s", normal_arr[i];
    print "";
}
{
    event_id = $1;
    split($2, psi1, ",");
    split($3, psi2, ",");
    printf "%s", event_id;
    for(i=1;i<=ntumor;i++) printf "\t%s", (i in psi1)? psi1[i] : "NA";
    for(i=1;i<=nnormal;i++) printf "\t%s", (i in psi2)? psi2[i] : "NA";
    print "";
}
' "$input_file" > "$output_file"

echo "✅ PSI matrix annotated and saved to: $output_file"