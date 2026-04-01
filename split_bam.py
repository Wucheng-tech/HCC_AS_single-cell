#!/usr/bin/env python3
import pysam
import os
import gzip
import pandas as pd
from multiprocessing import Pool
import argparse

def load_barcodes(barcode_file):
    """
    加载有效的 cell barcodes：
    - 如果是 barcodes.tsv.gz，不跳过；
    - 否则跳过第一行（比如 .xls 带标题的）。
    """
    basename = os.path.basename(barcode_file)
    skip_header = 0 if basename.endswith("barcodes.tsv.gz") else 1
    barcodes = pd.read_csv(barcode_file, sep="\t", skiprows=skip_header, header=None, usecols=[0])
    return barcodes.iloc[:, 0].tolist()

def process_barcode(args):
    """
    子进程函数：根据 barcode 过滤 BAM 文件中的 reads 写入对应文件
    """
    barcode, input_bam, output_dir = args
    output_bam_path = os.path.join(output_dir, f"{barcode}.bam")
    with pysam.AlignmentFile(input_bam, "rb") as bam_file, \
         pysam.AlignmentFile(output_bam_path, "wb", header=bam_file.header) as output_bam:
        for read in bam_file.fetch(until_eof=True):
            # 取 read name 的前缀作为匹配依据，常见 10x 格式为：BARCODE_SOMEINFO
            read_barcode = read.query_name.split("_")[0]
            if read_barcode == barcode:
                output_bam.write(read)

def main(sample_id, threads=10):
    """
    主函数，根据 sample ID 推导路径，读取 barcode 列表并启动多进程处理
    """
    base_path = "/public/home/wucheng/Analysis/scFAST/Data"
    bam_file = os.path.join(base_path, sample_id, "step2", "STAR", f"{sample_id}_SortedByCoordinate.bam")
    barcode_file = os.path.join(base_path, sample_id, "step3", "filtered_feature_bc_matrix", "barcodes.tsv.gz")
    output_dir = os.path.join(base_path, sample_id, "step2", "STAR", "split_by_barcode")

    os.makedirs(output_dir, exist_ok=True)

    print(f"📦 加载 barcodes from: {barcode_file}")
    barcodes = load_barcodes(barcode_file)
    print(f"✅ 共加载 {len(barcodes)} 个 barcodes")

    # 创建任务参数列表
    args_list = [(barcode, bam_file, output_dir) for barcode in barcodes]

    print(f"🚀 使用 {threads} 线程开始处理 BAM 文件...")
    with Pool(threads) as pool:
        pool.map(process_barcode, args_list)
    print("🎉 拆分完成！")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="根据 sample ID 拆分 BAM 文件为每个 barcode 的子文件")
    parser.add_argument("sample_id", help="样本 ID，例如 P243091_M")
    parser.add_argument("-t", "--threads", type=int, default=10, help="并行线程数（默认10）")
    args = parser.parse_args()

    main(args.sample_id, args.threads)
