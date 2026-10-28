#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pysam

def parse_args():
    p = argparse.ArgumentParser(
        description="统计 BAM 中一对 reads 分别位于 chr1 和 chr2 的数量"
    )
    p.add_argument("--bam", required=True,
                   help="输入 BAM 文件（已排序并建立索引）")
    p.add_argument("--chr1", required=True,
                   help="染色体1名称，例如 Chr1")
    p.add_argument("--chr2", required=True,
                   help="染色体2名称，例如 Chr2")
    p.add_argument("--proper_only", action="store_true",
                   help="只统计 proper pair（read.is_proper_pair 为 True）")
    return p.parse_args()

def main():
    args = parse_args()

    bam = pysam.AlignmentFile(args.bam, "rb")

    chr1_name = args.chr1
    chr2_name = args.chr2

    total_pairs = 0
    pairs_chr1_chr2 = 0
    pairs_chr1_chr2_dir1 = 0  # read1 在 chr1，read2 在 chr2
    pairs_chr1_chr2_dir2 = 0  # read1 在 chr2，read2 在 chr1

    for read in bam.fetch(until_eof=True):
        # 只统计成对 reads 的 read1，避免重复计数
        if not read.is_paired or not read.is_read1:
            continue
        if read.is_unmapped or read.mate_is_unmapped:
            continue
        if args.proper_only and (not read.is_proper_pair):
            continue

        rchr = bam.get_reference_name(read.reference_id)
        mchr = bam.get_reference_name(read.next_reference_id)

        total_pairs += 1

        # 一端在 chr1，一端在 chr2
        if ((rchr == chr1_name and mchr == chr2_name) or
            (rchr == chr2_name and mchr == chr1_name)):
            pairs_chr1_chr2 += 1
            if rchr == chr1_name and mchr == chr2_name:
                pairs_chr1_chr2_dir1 += 1
            else:
                pairs_chr1_chr2_dir2 += 1

    bam.close()

    print(f"# BAM: {args.bam}")
    print(f"# chr1: {chr1_name}, chr2: {chr2_name}")
    print(f"total_read_pairs\t{total_pairs}")
    print(f"pairs_chr1_chr2_total\t{pairs_chr1_chr2}")
    print(f"pairs_chr1(chr1)-chr2(chr2)\t{pairs_chr1_chr2_dir1}")
    print(f"pairs_chr2(chr2)-chr1(chr1)\t{pairs_chr1_chr2_dir2}")

if __name__ == "__main__":
    main()
