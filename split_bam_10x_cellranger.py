#!/usr/bin/env python
# -*- coding: utf-8 -*-
#source activate /public/home/jiangzh/miniconda3/envs/scape_2
#python split_bam_10x_cellranger.py ./outs -t 20
import pysam
import os
import pandas as pd
from multiprocessing import Pool
import argparse


def load_barcodes(barcode_file):
    """
    Read 10x barcodes.tsv.gz file
    Each line contains one cell barcode, e.g.:
    AAACCCAAGAAACCAT-1
    """
    barcodes = pd.read_csv(barcode_file, sep="\t", header=None)

    return barcodes.iloc[:, 0].tolist()


def process_barcode(args):
    """
    Split BAM file by CB tag (cell barcode)

    For each barcode:
    extract reads with matching CB tag
    and write to a separate BAM file
    """

    barcode, input_bam, output_dir = args

    output_bam = os.path.join(
        output_dir,
        "{}.bam".format(barcode)   # compatible with old python
    )

    with pysam.AlignmentFile(input_bam, "rb") as bam, \
         pysam.AlignmentFile(output_bam, "wb", header=bam.header) as out_bam:

        for read in bam.fetch(until_eof=True):

            if read.has_tag("CB"):

                if read.get_tag("CB") == barcode:

                    out_bam.write(read)


def main(outs_path, threads):

    # path to 10x BAM file
    bam_file = os.path.join(
        outs_path,
        "possorted_genome_bam.bam"
    )

    # path to barcode list
    barcode_file = os.path.join(
        outs_path,
        "filtered_feature_bc_matrix",
        "barcodes.tsv.gz"
    )

    # output directory
    output_dir = os.path.join(
        outs_path,
        "split_by_barcode"
    )

    os.makedirs(output_dir, exist_ok=True)

    print("outs path:", outs_path)
    print("BAM file:", bam_file)
    print("barcode file:", barcode_file)

    barcodes = load_barcodes(barcode_file)

    print("Total cells detected: {}".format(len(barcodes)))

    args_list = [

        (bc, bam_file, output_dir)

        for bc in barcodes
    ]

    print("Start splitting BAM (threads={})".format(threads))

    with Pool(threads) as pool:

        pool.map(process_barcode, args_list)

    print("Finished")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Split 10x BAM into one BAM per cell using CB tag"
    )

    parser.add_argument(
        "outs_path",
        help="Path to 10x outs directory, e.g. /path/to/outs"
    )

    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=10,
        help="Number of threads (default=10)"
    )

    args = parser.parse_args()

    main(args.outs_path, args.threads)