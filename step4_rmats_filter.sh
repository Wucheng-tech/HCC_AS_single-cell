#!/bin/bash

############################################################
# Step4: Post-process rmats results with filtering
# Usage: bash step4_rmats_filter.sh /path/to/parent_dir
# 1. The script will check all subdirectories under the parent directory
# 2. If summary.txt is complete -> run rmats_filtering.py
# 3. If summary.txt is incomplete -> rerun rmats_jc_postprocess.sh to generate .new.txt
# 4. Then run rmats_filtering.py
############################################################

# Check input argument
if [[ -z "$1" ]]; then
    echo "Usage: bash $0 /path/to/parent_dir"
    exit 1
fi

parent_dir=$1

# Activate rmats environment
source activate /public/home/jiangzh/miniconda3/envs/rmats

# Paths
rmats_filter=/public/home/wucheng/Analysis/AS/bulk_seq/software/rmats-turbo-tutorial-main/scripts/rmats_filtering.py
rmats_postprocess=/public/home/wucheng/Analysis/AS/scRNA/scFAST/Sample/AA/rmats_jc_postprocess.sh

# Event types
event_array=("SE" "MXE" "RI" "A5SS" "A3SS")

# Find all subdirectories named 'Tumor_Normal' under the parent directory
mapfile -t base_dirs < <(find "$parent_dir" -type d -path "*/merged/T_N/Tumor_Normal")

if [[ ${#base_dirs[@]} -eq 0 ]]; then
    echo "❌ No directories found under $parent_dir"
    exit 1
fi

# Iterate over directories
for dir in "${base_dirs[@]}"; do
    echo "📂 Processing directory: $dir"
    
    # Ensure Fit_new exists
    mkdir -p "${dir}/Fit_new"
    cd "${dir}/Fit_new" || { echo "❌ Failed to cd $dir/Fit_new"; continue; }

    # Check summary.txt for completeness
    summary_file="${dir}/summary.txt"
    complete_flag=true
    if [[ ! -f "$summary_file" ]]; then
        echo "⚠️ summary.txt not found, will regenerate .new.txt"
        complete_flag=false
    else
        # Optional: check number of lines or expected content
        line_count=$(wc -l < "$summary_file")
        if (( line_count < 5 )); then
            echo "⚠️ summary.txt incomplete (lines=$line_count), will regenerate .new.txt"
            complete_flag=false
        fi
    fi

    # If incomplete -> rerun rmats postprocess
    if [[ "$complete_flag" == false ]]; then
        echo "🔄 Running rmats_jc_postprocess.sh to generate .new.txt..."
        bash "$rmats_postprocess" "$dir"
    fi

    # Run rmats_filtering.py for all events
    echo "📋 Running rmats_filtering.py for events..."
    for event in "${event_array[@]}"; do
        input_file=""
        if [[ -f "${dir}/${event}.MATS.JC.new.txt" ]]; then
            input_file="${dir}/${event}.MATS.JC.new.txt"
        elif [[ -f "${dir}/${event}.MATS.JC.txt" ]]; then
            input_file="${dir}/${event}.MATS.JC.txt"
        else
            echo "⚠️ Event file not found for $event, skipping..."
            continue
        fi
        
        python "$rmats_filter" "$input_file" 5,0.05,0.95,0.01,0.05,0.5,0.2
    done

    echo "✅ Finished processing $dir"
done

echo "🎉 Step4: rmats postprocessing and filtering completed."