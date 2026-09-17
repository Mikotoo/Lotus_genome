#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pysam


def in_region(chrom, pos1based, target_chr, start, end):
    """判断 (chrom, pos) 是否落在 target_chr:start-end (1-based, 闭区间)"""
    return (chrom == target_chr) and (start <= pos1based <= end)


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "在参考基因组坐标下统计：一端在 Chr1[1-bp]，"
            "另一端在 Chr1[bp+1-end] 或 Chr2 任意位置的成对 reads 数量"
        )
    )
    p.add_argument("--bam", required=True,
                   help="比对到 Gifu 基因组的 BAM（已 sort + index）")
    p.add_argument("--chr1", default="Chr1",
                   help="Chr1 名称（默认 Chr1）")
    p.add_argument("--chr2", default="Chr2",
                   help="Chr2 名称（默认 Chr2）")
    p.add_argument("--bp", type=int, required=True,
                   help="Chr1 上的断点坐标（例如 22839936）")
    p.add_argument("--min-mapq", type=int, default=0,
                   help="最小比对质量（MAPQ），默认 0 不过滤")
    return p.parse_args()


def main():
    args = parse_args()
    bam = pysam.AlignmentFile(args.bam, "rb")

    chr1 = args.chr1
    chr2 = args.chr2
    bp   = args.bp

    # 1) 检查染色体名是否存在
    try:
        chr1_len = bam.get_reference_length(chr1)
    except ValueError:
        raise SystemExit(f"[ERROR] 在 BAM header 中找不到 {chr1}，"
                         f"实际染色体名可能是：{bam.references[:10]} ...")
    try:
        chr2_len = bam.get_reference_length(chr2)
    except ValueError:
        raise SystemExit(f"[ERROR] 在 BAM header 中找不到 {chr2}，"
                         f"实际染色体名可能是：{bam.references[:10]} ...")

    # 2) 定义区域（1-based 闭区间）
    region1_chr   = chr1
    region1_start = 1
    region1_end   = bp

    region2a_chr   = chr1
    region2a_start = bp + 1
    region2a_end   = chr1_len

    region2b_chr   = chr2
    region2b_start = 1
    region2b_end   = chr2_len

    total_pairs = 0
    n_1_to_1tail = 0   # 区域1 <-> Chr1 tail
    n_1_to_chr2  = 0   # 区域1 <-> Chr2
    n_read1_in_region1 = 0  # 有多少 read1 落在 region1（自检用）

    # 只看 read1，避免重复计数
    for read in bam.fetch(until_eof=True):
        if not read.is_paired or not read.is_read1:
            continue
        if read.is_unmapped or read.mate_is_unmapped:
            continue
        if read.is_secondary or read.is_supplementary:
            continue
        if read.mapping_quality < args.min_mapq:
            continue

        rchr = bam.get_reference_name(read.reference_id)
        rpos = read.reference_start + 1  # 1-based

        mate_id = read.next_reference_id
        if mate_id < 0:
            continue
        mchr = bam.get_reference_name(mate_id)
        mpos = read.next_reference_start + 1  # 1-based

        total_pairs += 1

        r_in_1     = in_region(rchr, rpos, region1_chr,   region1_start,  region1_end)
        r_in_1tail = in_region(rchr, rpos, region2a_chr,  region2a_start, region2a_end)
        r_in_chr2  = in_region(rchr, rpos, region2b_chr,  region2b_start, region2b_end)

        m_in_1     = in_region(mchr, mpos, region1_chr,   region1_start,  region1_end)
        m_in_1tail = in_region(mchr, mpos, region2a_chr,  region2a_start, region2a_end)
        m_in_chr2  = in_region(mchr, mpos, region2b_chr,  region2b_start, region2b_end)

        if r_in_1:
            n_read1_in_region1 += 1

        # 一端在区域1，另一端在 Chr1 tail
        if (r_in_1 and m_in_1tail) or (r_in_1tail and m_in_1):
            n_1_to_1tail += 1

        # 一端在区域1，另一端在 Chr2
        if (r_in_1 and m_in_chr2) or (r_in_chr2 and m_in_1):
            n_1_to_chr2 += 1

    bam.close()

    print(f"# BAM: {args.bam}")
    print(f"# Chr1: {chr1} (len={chr1_len}), Chr2: {chr2} (len={chr2_len})")
    print(f"# Chr1 break at: {bp}")
    print(f"# region1  = {chr1}:{region1_start}-{region1_end}")
    print(f"# region2a = {chr1}:{region2a_start}-{region2a_end} (Chr1 tail)")
    print(f"# region2b = {chr2}:{region2b_start}-{region2b_end} (Chr2 full)")
    print(f"total_read_pairs\t{total_pairs}")
    print(f"read1_in_region1\t{n_read1_in_region1}")
    print(f"pairs_region1_to_chr1tail\t{n_1_to_1tail}")
    print(f"pairs_region1_to_chr2\t{n_1_to_chr2}")


if __name__ == "__main__":
    main()
