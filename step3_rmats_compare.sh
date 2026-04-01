#!/bin/bash

############################################################
# Usage:
# bash step3_rmats_compare.sh /path/to/merged_bams CompareGroupName
# Example:
# bash step3_rmats_compare.sh /public/home/wucheng/Analysis/AS/scRNA/scFAST/Sample/AA/celltype_bin_100/merged Tumor_Normal
############################################################

# Input arguments
input_dir="$1"        # Directory containing merged BAM files
compare_group="$2"    # Comparison group name, e.g., Tumor_Normal
read_length=150       # Read length for rmats
threads=20            # Number of threads for rmats
gtf_file="/public/home/wucheng/software/RNAediting/reference/gencode.v44.primary_assembly.annotation.gtf"
rmats_bin="/public/home/jiangzh/miniconda3/envs/rmats/bin/rmats.py"

# Check input
if [[ -z "$input_dir" || -z "$compare_group" ]]; then
    echo "Usage:"
    echo "bash step3_rmats_compare.sh /path/to/merged_bams CompareGroupName"
    exit 1
fi

# Output directories
output_dir="${input_dir}/T_N"
mkdir -p "$output_dir"

tumor_out="${output_dir}/Tumor_BAM_fls.txt"
normal_out="${output_dir}/Normal_BAM_fls.txt"

tmp_tumor="${tumor_out}.tmp"
tmp_normal="${normal_out}.tmp"

# Clear temporary files if they exist
> "$tmp_tumor"
> "$tmp_normal"

############################################################
# Collect BAM files into tumor and normal lists
############################################################
echo "📋 Generating BAM lists..."
find "$input_dir" -maxdepth 1 -name "*.bam" | while read -r bam_file; do
    filename=$(basename "$bam_file")
    if [[ "$filename" == *_PN_* ]]; then
        # Normal sample
        echo "$bam_file" >> "$tmp_normal"
    elif [[ "$filename" == *_T1_* || "$filename" == *_T2_* ]]; then
        # Tumor sample
        echo "$bam_file" >> "$tmp_tumor"
    fi
done

# Sort and join BAM paths with commas (single line, no trailing comma)
tumor_line=$(sort "$tmp_tumor" | paste -sd, -)
normal_line=$(sort "$tmp_normal" | paste -sd, -)
echo "$tumor_line" > "$tumor_out"
echo "$normal_line" > "$normal_out"

# Remove temporary files
rm "$tmp_tumor" "$tmp_normal"

echo "✅ BAM lists created:"
echo "  Tumor: $tumor_out"
echo "  Normal: $normal_out"

############################################################
# Run rmats analysis
############################################################
# Activate rmats conda environment
source activate /public/home/jiangzh/miniconda3/envs/rmats

rmats_out="${output_dir}/${compare_group}"
rmats_tmp="${rmats_out}_tmp_rmats"
mkdir -p "$rmats_tmp"

echo "🔄 Running rmats.py for $compare_group ..."
python "$rmats_bin" \
    --b1 "$tumor_out" \
    --b2 "$normal_out" \
    --gtf "$gtf_file" \
    --od "$rmats_out" \
    --tmp "$rmats_tmp" \
    -t single \
    --variable-read-length \
    --readLength "$read_length" \
    --nthread "$threads"

echo "🎉 rmats analysis completed. Results in: $rmats_out"
