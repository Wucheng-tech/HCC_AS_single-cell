#!/bin/bash

############################################################
# Step5: Generate PSI matrix from Fit_new directory
# Input: Fit_new directory containing filtered rmats files
# Output: IncLevel_summary.txt, IncLevel_summary_annotated.txt, MetaInfo_summary.txt
# Tumor/Normal BAM list inferred automatically from the path
############################################################

fit_new_dir="$1"  # Fit_new directory
# Example: /public/home/wucheng/Analysis/AS/scRNA/scFAST/Sample/AA/celltype_split_100/merged/T_N/Tumor_Normal/Fit_new

if [[ -z "$fit_new_dir" ]]; then
    echo "Usage: bash $0 /path/to/Fit_new"
    exit 1
fi

echo "📋 Processing Fit_new directory: $fit_new_dir"

# Infer parent directories
parent_dir=$(dirname "$fit_new_dir")  # .../Tumor_Normal
tn_dir=$(dirname "$parent_dir")       # .../T_N
tumor_file="$tn_dir/Tumor_BAM_fls.txt"
normal_file="$tn_dir/Normal_BAM_fls.txt"

if [[ ! -f "$tumor_file" || ! -f "$normal_file" ]]; then
    echo "❌ Tumor or Normal BAM list not found:"
    echo "Tumor: $tumor_file"
    echo "Normal: $normal_file"
    exit 1
fi

echo "✅ Found BAM lists:"
echo "Tumor: $tumor_file"
echo "Normal: $normal_file"

# Output files
inclevel_file="$fit_new_dir/IncLevel_summary.txt"
annotated_file="$fit_new_dir/IncLevel_summary_annotated.txt"
meta_file="$fit_new_dir/MetaInfo_summary.txt"

> "$inclevel_file"
> "$annotated_file"
> "$meta_file"

event_array=("SE" "RI" "MXE" "A5SS" "A3SS")

############################################################
# 1. Combine filtered rmats results to construct IncLevel summary
############################################################
for event in "${event_array[@]}"; do
    # Prefer .new.txt if exists
    if [[ -f "$fit_new_dir/filtered_${event}.MATS.JC.new.txt" ]]; then
        input_file="$fit_new_dir/filtered_${event}.MATS.JC.new.txt"
    elif [[ -f "$fit_new_dir/filtered_${event}.MATS.JC.txt" ]]; then
        input_file="$fit_new_dir/filtered_${event}.MATS.JC.txt"
    else
        echo "⚠️  File not found for event $event, skipping..."
        continue
    fi

    echo "✅ Processing: $input_file"

    # Get column indices for IncLevel1 and IncLevel2
    header=$(head -n 1 "$input_file")
    IFS=$'\t' read -r -a cols <<< "$header"
    inc1_idx=0
    inc2_idx=0
    for i in "${!cols[@]}"; do
        [[ "${cols[$i]}" == "IncLevel1" ]] && inc1_idx=$((i+1))
        [[ "${cols[$i]}" == "IncLevel2" ]] && inc2_idx=$((i+1))
    done

    if (( inc1_idx == 0 || inc2_idx == 0 )); then
        echo "⚠️ IncLevel columns not found for $event, skipping..."
        continue
    fi

    # Extract IncLevel1/2 columns with event identifier
    awk -v e="$event" -v i1="$inc1_idx" -v i2="$inc2_idx" 'BEGIN{FS=OFS="\t"; n=0}
    NR>1 {n++; print e"_"n, $i1, $i2}' "$input_file" >> "$inclevel_file"

    # ----------------------
    # Generate MetaInfo: event ID + first 5 columns
    awk -v e="$event" 'BEGIN{FS=OFS="\t"; n=0}
    NR>1 {n++; print e"_"n, $1, $2, $3, $4, $5}' "$input_file" >> "$meta_file"
done

echo "🎉 IncLevel summary written to: $inclevel_file"
echo "🎉 Meta info file created: $meta_file"

############################################################
# 2. Generate PSI matrix annotated with Tumor/Normal IDs
############################################################

# Parse Tumor BAM IDs
tumor_samples=()
IFS=',' read -ra tumor_paths <<< "$(cat "$tumor_file")"
for path in "${tumor_paths[@]}"; do
    id=$(basename "$path" .bam)
    tumor_samples+=("$id")
done

# Parse Normal BAM IDs
normal_samples=()
IFS=',' read -ra normal_paths <<< "$(cat "$normal_file")"
for path in "${normal_paths[@]}"; do
    id=$(basename "$path" .bam)
    normal_samples+=("$id")
done

# Write temporary files
printf "%s\n" "${tumor_samples[@]}" > "$fit_new_dir/tumor_ids.txt"
printf "%s\n" "${normal_samples[@]}" > "$fit_new_dir/normal_ids.txt"

ntumor=${#tumor_samples[@]}
nnormal=${#normal_samples[@]}

echo "📋 Generating annotated PSI matrix..."

awk -v ntumor="$ntumor" -v nnormal="$nnormal" '
BEGIN {
    FS=OFS="\t";
    # Load tumor IDs
    i=0; while ((getline line < "'$fit_new_dir'/tumor_ids.txt") > 0) {tumor_arr[++i]=line} close("'$fit_new_dir'/tumor_ids.txt");
    # Load normal IDs
    j=0; while ((getline line < "'$fit_new_dir'/normal_ids.txt") > 0) {normal_arr[++j]=line} close("'$fit_new_dir'/normal_ids.txt");
    # Print header
    printf "ID";
    for (k=1;k<=ntumor;k++) printf "\tPSI_Tumor_%s", tumor_arr[k];
    for (k=1;k<=nnormal;k++) printf "\tPSI_Normal_%s", normal_arr[k];
    print "";
}
{
    event_id=$1;
    split($2, psi1, ",");
    split($3, psi2, ",");
    printf "%s", event_id;
    for(i=1;i<=ntumor;i++) printf "\t%s", (i in psi1 ? psi1[i] : "NA");
    for(i=1;i<=nnormal;i++) printf "\t%s", (i in psi2 ? psi2[i] : "NA");
    print "";
}' "$inclevel_file" > "$annotated_file"

echo "🎉 PSI matrix generated: $annotated_file"