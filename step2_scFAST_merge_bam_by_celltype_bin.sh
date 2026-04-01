#!/bin/bash

############################################################
# Usage
# dos2unix step2_scFAST_merge_bam_by_celltype_bin.sh
# bash step2_scFAST_merge_bam_by_celltype_bin.sh 100
# bash step2_scFAST_merge_bam_by_celltype_bin.sh 150 random
# bash step2_scFAST_merge_bam_by_celltype_bin.sh all 
############################################################

bin_size=$1
mode=$2   # optional: random

if [[ -z "$bin_size" ]]; then
    echo "Usage:"
    echo "bash split_merge_bam_by_celltype_bin.sh 100"
    echo "bash split_merge_bam_by_celltype_bin.sh 100 random"
    echo "bash split_merge_bam_by_celltype_bin.sh all" ##celltype merge
    exit 1
fi

############################################################
# paths
############################################################

input_file="/public/home/wucheng/Analysis/AS/scRNA/scFAST/sample_cell.txt"
bam_base="/public/home/wucheng/Analysis/AS/scRNA/scFAST/Sample"
output_base="/public/home/wucheng/Analysis/AS/scRNA/scFAST/Sample/AA/celltype_bin_${bin_size}"
if [[ "$mode" == "random" ]]; then
    output_base="${output_base}_random"
fi
mkdir -p "$output_base"

############################################################
# collect BAM paths
############################################################

temp_file="/tmp/all_bam_paths_${bin_size}_${mode}.tsv"

> "$temp_file"
echo "📋 Collecting BAM paths..."
tail -n +2 "$input_file" | \
awk -F'"' '{print $2"\t"$4"\t"$8"\t"$10}' | \
sort | uniq | \
while IFS=$'\t' read -r Sample Celltype Cell ID
do
    bam_path="${bam_base}/${Sample}/${Cell}.bam"
    if [[ -f "$bam_path" ]]; then
        echo -e "${Sample}\t${Celltype}\t${bam_path}" >> "$temp_file"
    fi
done

############################################################
# split into bins or merge all
############################################################
cut -f1,2 "$temp_file" | sort | uniq |
while IFS=$'\t' read -r Sample Celltype
do
    # 所有该 Sample-Celltype 的 BAM 路径
    group_bams=$(awk -F'\t' -v s="$Sample" -v c="$Celltype" '$1==s && $2==c {print $3}' "$temp_file")
    count=$(echo "$group_bams" | wc -l)

    if [[ "$bin_size" == "all" ]]; then
        # Merge all BAMs
        merged_dir="${output_base}/merged"
        mkdir -p "$merged_dir"
        out_bam="${merged_dir}/${Sample}_${Celltype}.bam"
        echo "🔄 Merging ALL $count BAMs for $Sample-$Celltype → $out_bam"
        samtools merge -f "$out_bam" -b <(echo "$group_bams")
        samtools index "$out_bam"
    else
        if (( count >= bin_size )); then
            echo "✅ $Sample-$Celltype ($count cells) → splitting bins of $bin_size..."
            bams_to_use="$group_bams"
            [[ "$mode" == "random" ]] && bams_to_use=$(echo "$group_bams" | shuf)

            # split into bin_size lines
            echo "$bams_to_use" | split -l "$bin_size" - "$output_base/${Sample}_${Celltype}_"

            idx=1
            for chunk in "$output_base/${Sample}_${Celltype}_"*; do
                final_list="${output_base}/${Sample}_${Celltype}_${idx}.list"
                mv "$chunk" "$final_list"

                # 只有当.list行数==bin_size才进行合并
                lines_in_list=$(wc -l < "$final_list")
                if (( lines_in_list == bin_size )); then
                    merged_dir="${output_base}/merged"
                    mkdir -p "$merged_dir"
                    out_bam="${merged_dir}/$(basename "$final_list" .list).bam"
                    echo "🔄 Merging bin $(basename "$final_list") → $out_bam"
                    samtools merge -f "$out_bam" -b "$final_list"
                    samtools index "$out_bam"
                else
                    echo "⚠️ Bin $(basename "$final_list") has $lines_in_list cells < bin_size ($bin_size), skipping BAM merge"
                fi

                ((idx++))
            done
        else
            echo "⚠️ $Sample-$Celltype only $count cells, skipped."
        fi
    fi
done

echo "🎉 Done!"
echo "output directory: ${output_base}/merged"

############################################################
# Generate cell merge summary from .list files
############################################################

list_dir="$output_base"  # use same bin output directory
merge_summary="${output_base}/Cell_merge_summary.txt"
> "$merge_summary"  # clear old file if exists

echo "📋 Generating cell merge summary from .list files in $list_dir..."

for list_file in "$list_dir"/*.list; do
    base_name=$(basename "$list_file" .list)
    line_count=$(wc -l < "$list_file")
    # Only consider full bins matching bin_size (optional)
    if [[ "$bin_size" != "all" && "$line_count" -eq "$bin_size" ]] || [[ "$bin_size" == "all" ]]; then
        while read -r bam_path; do
            bam_name=$(basename "$bam_path" .bam)
            echo -e "${base_name}\t${bam_name}" >> "$merge_summary"
        done < "$list_file"
    fi
done

echo "✅ Cell merge summary written to: $merge_summary" 