#!/bin/bash

## =======================================================
## Huh7_SR vs Huh7 cell lines: rmats analysis & PSI matrix construction
## =======================================================

# -------------------------------
# 1. Define output BAM list files
# -------------------------------
Huh7_SR_out="/public/home/wucheng/Analysis/AS/Huh7/Huh7_SR_BAM_list.txt"
Huh7_out="/public/home/wucheng/Analysis/AS/Huh7/Huh7_BAM_list.txt"

# Sample IDs for Huh7 and Huh7_SR
samples=("Huh7_1" "Huh7_2" "Huh7_3" "Huh7_SR_1" "Huh7_SR_2" "Huh7_SR_3")

# -------------------------------
# 2. Generate BAM file lists
# -------------------------------
for i in "${samples[@]}"
do
    bam_path="/public/home/wucheng/Analysis/AS/Huh7/Bam/${i}.bam"

    if [[ "$i" == Huh7_SR_* ]]; then
        echo -n "$bam_path," >> "$Huh7_SR_out"
    elif [[ "$i" == Huh7_[0-9]* ]]; then
        echo -n "$bam_path," >> "$Huh7_out"
    fi
done

# Remove trailing commas
sed -i 's/,\s*$//' "$Huh7_SR_out"
sed -i 's/,\s*$//' "$Huh7_out"

# -------------------------------
# 3. Activate rmats conda environment
# -------------------------------
source activate /public/home/jiangzh/miniconda3/envs/rmats

# -------------------------------
# 4. Run rmats for Huh7_SR vs Huh7
# -------------------------------
compare_group="Huh7_SR_Huh7"
/public/home/jiangzh/miniconda3/envs/rmats/bin/rmats.py \
    --b1 "$Huh7_SR_out" \
    --b2 "$Huh7_out" \
    --gtf /public/home/wucheng/software/RNAediting/reference/gencode.v44.primary_assembly.annotation.gtf \
    --od /public/home/wucheng/Analysis/AS/Huh7/${compare_group} \
    --tmp /public/home/wucheng/Analysis/AS/Huh7/${compare_group}_tmp_rmats \
    -t paired \
    --variable-read-length \
    --readLength 150 \
    --nthread 20

# =======================================================
# 5. Filtering rmats results & differential splicing
# =======================================================
rmats_filter=/public/home/wucheng/Analysis/AS/bulk_seq/software/rmats-turbo-tutorial-main/scripts/rmats_filtering.py
declare -a event_array=("SE" "MXE" "RI" "A5SS" "A3SS")

# Enter the rmats output directory
cd /public/home/wucheng/Analysis/AS/Huh7/${compare_group}/Fit

# Filter each AS event type
for event in "${event_array[@]}"; do
    python $rmats_filter /public/home/wucheng/Analysis/AS/Huh7/${compare_group}/${event}.MATS.JC.txt 10,0.05,0.95,0.01,0.05,0.5,0.2
done

# =======================================================
# 6. Merge all events into an RSI (Inclusion Level) matrix
# =======================================================
input_dir="/public/home/wucheng/Analysis/AS/Huh7/${compare_group}/Fit"
output_file="${input_dir}/IncLevel_summary.txt"

# Clear output file
> "$output_file"

# Extract IncLevel1 and IncLevel2 for each filtered event
for event in SE RI MXE A3SS A5SS; do
    input_file="${input_dir}/filtered_${event}.MATS.JC.txt"
    if [[ ! -f "$input_file" ]]; then
        echo "❌ File not found: $input_file"
        continue
    fi
    echo "✅ Processing: $input_file"

    # Find column indices for IncLevel1 and IncLevel2
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

    # Extract and append with event ID
    awk -v e="$event" -v i1="$inc1_idx" -v i2="$inc2_idx" 'BEGIN{FS=OFS="\t"; n=0}
    NR>1 { n++; print e"_"n, $i1, $i2 }' "$input_file" >> "$output_file"
done

echo "🎉 Done! Combined IncLevel file written to: $output_file"

# =======================================================
# 7. Construct meta-information file
# =======================================================
output_meta="${input_dir}/MetaInfo_summary.txt"
> "$output_meta" # Clear output file

for event in SE RI MXE A3SS A5SS; do
    input_file="${input_dir}/filtered_${event}.MATS.JC.txt"
    if [[ ! -f "$input_file" ]]; then
        echo "❌ File not found: $input_file"
        continue
    fi
    echo "✅ Processing: $input_file"

    # Extract event ID and first five columns
    awk -v e="$event" 'BEGIN{FS=OFS="\t"; n=0}
    NR==1 { next }
    { n++; print e"_"n, $1, $2, $3, $4, $5 }' "$input_file" >> "$output_meta"
done

echo "🎉 Meta-information file created: $output_meta"

# =======================================================
# 8. Generate PSI expression matrix (annotated)
# =======================================================
input_file="${input_dir}/IncLevel_summary.txt"
output_file="${input_dir}/IncLevel_summary_annotated.txt"
tumor_file="$Huh7_SR_out"
normal_file="$Huh7_out"

# Parse tumor sample IDs
tumor_samples=()
IFS=',' read -ra tumor_paths <<< "$(cat "$tumor_file")"
for path in "${tumor_paths[@]}"; do
    rel_path="${path##*/celltype/}"
    id=$(echo "$rel_path" | sed 's#/#__#g')
    tumor_samples+=("$id")
done

# Parse normal sample IDs
normal_samples=()
IFS=',' read -ra normal_paths <<< "$(cat "$normal_file")"
for path in "${normal_paths[@]}"; do
    rel_path="${path##*/celltype/}"
    id=$(echo "$rel_path" | sed 's#/#__#g')
    normal_samples+=("$id")
done

# Convert arrays to comma-separated strings
ntumor=${#tumor_samples[@]}
nnormal=${#normal_samples[@]}
tumor_str=$(IFS=,; echo "${tumor_samples[*]}")
normal_str=$(IFS=,; echo "${normal_samples[*]}")

# Annotate PSI matrix with sample IDs
awk -v tumors="$tumor_str" -v normals="$normal_str" -v ntumor="$ntumor" -v nnormal="$nnormal" '
BEGIN {
    FS=OFS="\t";
    split(tumors, tumor_arr, ",");
    split(normals, normal_arr, ",");
    # Print header
    printf "ID";
    for (i=1; i<=ntumor; i++) printf "\tPSI_Tumor_%s", tumor_arr[i];
    for (i=1; i<=nnormal; i++) printf "\tPSI_Normal_%s", normal_arr[i];
    print "";
}
{
    event_id = $1;
    split($2, psi1, ",");
    split($3, psi2, ",");
    printf "%s", event_id;
    for (i=1; i<=ntumor; i++) printf "\t%s", (i in psi1)? psi1[i] : "NA";
    for (i=1; i<=nnormal; i++) printf "\t%s", (i in psi2)? psi2[i] : "NA";
    print "";
}
' "$input_file" > "$output_file"

echo "✅ PSI expression matrix annotated and saved to: $output_file"