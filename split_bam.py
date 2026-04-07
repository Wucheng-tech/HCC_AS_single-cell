#!/usr/bin/env python3
import pysam
import os
import gzip
import pandas as pd
from multiprocessing import Pool
import argparse

def load_barcodes(barcode_file):
    """
   Load valid cell barcodes.
   Behavior:
    - If the input file is 'barcodes.tsv.gz' (10x format), do not skip header.
    - Otherwise (e.g. Excel-exported .xls with header), skip the first line.
    """
    basename = os.path.basename(barcode_file)
    skip_header = 0 if basename.endswith("barcodes.tsv.gz") else 1
    barcodes = pd.read_csv(barcode_file, sep="\t", skiprows=skip_header, header=None, usecols=[0])
    return barcodes.iloc[:, 0].tolist()

def process_barcode(args):
    """
    Worker function executed in parallel.
    """
    barcode, input_bam, output_dir = args
    output_bam_path = os.path.join(output_dir, f"{barcode}.bam")
    with pysam.AlignmentFile(input_bam, "rb") as bam_file, \
         pysam.AlignmentFile(output_bam_path, "wb", header=bam_file.header) as output_bam:
        for read in bam_file.fetch(until_eof=True):
            # Extract barcode prefix from read name
            read_barcode = read.query_name.split("_")[0]
            if read_barcode == barcode:
                output_bam.write(read)

def main(sample_id, threads=10):
    """
    Main function.
    """
    base_path = "/public/home/wucheng/Analysis/scFAST/Data"
    bam_file = os.path.join(base_path, sample_id, "step2", "STAR", f"{sample_id}_SortedByCoordinate.bam")
    barcode_file = os.path.join(base_path, sample_id, "step3", "filtered_feature_bc_matrix", "barcodes.tsv.gz")
    output_dir = os.path.join(base_path, sample_id, "step2", "STAR", "split_by_barcode")

    os.makedirs(output_dir, exist_ok=True)

    print(f"Loading barcodes from: {barcode_file}")
    barcodes = load_barcodes(barcode_file)
    print(f"Total barcodes loaded: {len(barcodes)}")

    # create argument list for parallel processing
    args_list = [(barcode, bam_file, output_dir) for barcode in barcodes]

    print(f"Start splitting BAM file using {threads} processes...")
    with Pool(threads) as pool:
        pool.map(process_barcode, args_list)
    print("BAM splitting completed.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split BAM file into barcode-specific BAM files based on sample ID")
    parser.add_argument("sample_id", help="Sample ID (example: P243091_M)")
    parser.add_argument("-t", "--threads", type=int, default=10, help="Number of parallel processes (default: 10)")
    args = parser.parse_args()

    main(args.sample_id, args.threads)
