#!/bin/bash

############################################################
# Step 1
# Split STAR BAM into single-cell BAM
############################################################
base_dir="/public/home/wucheng/Analysis/scFAST/Data"
threads=20
samples=(
P267835_T2
P267835_T1
P267835_N
P267661_T1
P267661_T2
P267661_N
P264663_T
P264663_N
P264663_M
P262103_T1
P262103_N
P262103_M1
P243091_T
P243091_PN
P243091_M
)
for s in "${samples[@]}"
do
    work_dir="${base_dir}/${s}/step2/STAR"
    log_file="${work_dir}/split_${s}.log"
    echo "Splitting BAM for $s"
    cd "$work_dir"
    nohup python ${base_dir}/split_bam.py \
        "$s" \
        -t $threads \
        > "$log_file" 2>&1 &
done
echo "All split jobs submitted"

# organize cell BAM by Sample 
############################################################ 
input_file="/public/home/wucheng/Analysis/AS/scRNA/scFAST/Sample/sample_cell.txt" 
bam_base="/public/home/wucheng/Analysis/scFAST/Data" 
target_base="/public/home/wucheng/Analysis/AS/scRNA/scFAST/Sample" 
tail -n +2 "$input_file" | awk -F'"' '{print $2, $8, $10}' | while read Sample cell ID; do
    bam_file="${bam_base}/${ID}/step2/STAR/split_by_barcode/${cell}.bam"
    target_dir="${target_base}/${Sample}"
    if [ -f "$bam_file" ]; then
        mkdir -p "$target_dir"
        cp "$bam_file" "$target_dir/"
        echo "Copied: $bam_file -> $target_dir/"
    else
        echo "Missing: $bam_file" >&2
    fi
done
echo "✅ All done!"
