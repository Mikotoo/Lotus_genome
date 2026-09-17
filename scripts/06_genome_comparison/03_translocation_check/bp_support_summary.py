#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pysam

def parse_bp(s):
    """
    解析形如: Chr1:22839936:-  或 Chr2:90335659:+
    返回: (chrom, pos(int), side)
    side 为 '+' 或 '-'
    """
    parts = s.split(":")
    if len(parts) != 3:
        raise ValueError("断点格式必须为 Chr:Pos:+/- ，例如 Chr1:22839936:-")
    chrom = parts[0]
    pos = int(parts[1])
    side = parts[2]
    if side not in ["+", "-"]:
        raise ValueError("方向必须是 + 或 -")
    return chrom, pos, side

def bp_to_region(chrom, pos, side, win):
    """
    根据断点和方向生成 1-based 闭区间 [start, end]
    side == '-'  -> 左侧窗口: [pos - (win-1), pos]
    side == '+'  -> 右侧窗口: [pos, pos + (win-1)]
    """
    if side == "-":
        start = max(1, pos - (win - 1))
        end = pos
    else:  # '+'
        start = pos
        end = pos + (win - 1)
    return chrom, start, end

def in_region(chrom, pos1based, target_chr, start, end):
    return (chrom == target_chr) and (start <= pos1based <= end)

def count_pairs_between_regions(bam, chr1, s1, e1, chr2, s2, e2):
    """
    统计：一端在 region1，一端在 region2 的成对 reads 数（按 read name 去重）
    region1 = chr1:s1-e1, region2 = chr2:s2-e2
    """
    pair_names = set()

    def scan_region(this_chr, this_start, this_end,
                    other_chr, other_start, other_end):
        start0 = this_start - 1  # 1-based -> 0-based
        end0 = this_end          # pysam.fetch: [start0, end0)
        for read in bam.fetch(this_chr, start0, end0):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if not read.is_paired or read.mate_is_unmapped:
                continue

            pos1based = read.reference_start + 1
            if not (this_start <= pos1based <= this_end):
                continue

            mate_id = read.next_reference_id
            if mate_id < 0:
                continue
            mate_chr = bam.get_reference_name(mate_id)
            mate_pos = read.next_reference_start + 1

            if in_region(mate_chr, mate_pos, other_chr, other_start, other_end):
                pair_names.add(read.query_name)

    scan_region(chr1, s1, e1, chr2, s2, e2)
    scan_region(chr2, s2, e2, chr1, s1, e1)

    return len(pair_names)

def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "同时统计：一套 reads 在一个 BAM 中对真结构/假结构的支持 "
            "(true & false breakpoint pairs)."
        )
    )
    p.add_argument("--bam", required=True,
                   help="输入 BAM（排序并建立索引）")

    # 真结构的两个断点
    p.add_argument("--true-bp1", required=True,
                   help="真结构断点1: Chr:Pos:+/-  例如 Chr1:22839936:-")
    p.add_argument("--true-bp2", required=True,
                   help="真结构断点2: Chr:Pos:+/-")

    # 假结构的两个断点
    p.add_argument("--false-bp1", required=True,
                   help="假结构断点1: Chr:Pos:+/-")
    p.add_argument("--false-bp2", required=True,
                   help="假结构断点2: Chr:Pos:+/-")

    p.add_argument("--win", type=int, default=10000,
                   help="窗口大小（bp），默认 10000 = 10 kb")
    p.add_argument("--no-header", action="store_true",
                   help="不输出表头（用于 >> 追加到同一个 summary）")
    p.add_argument("--out", default="-",
                   help="输出文件名，默认 - 为 stdout")
    return p.parse_args()

def main():
    args = parse_args()
    bam = pysam.AlignmentFile(args.bam, "rb")

    # 解析真结构断点
    t_chr1_bp, t_pos1_bp, t_side1 = parse_bp(args.true_bp1)
    t_chr2_bp, t_pos2_bp, t_side2 = parse_bp(args.true_bp2)
    t_chr1, t_s1, t_e1 = bp_to_region(t_chr1_bp, t_pos1_bp, t_side1, args.win)
    t_chr2, t_s2, t_e2 = bp_to_region(t_chr2_bp, t_pos2_bp, t_side2, args.win)

    # 解析假结构断点
    f_chr1_bp, f_pos1_bp, f_side1 = parse_bp(args.false_bp1)
    f_chr2_bp, f_pos2_bp, f_side2 = parse_bp(args.false_bp2)
    f_chr1, f_s1, f_e1 = bp_to_region(f_chr1_bp, f_pos1_bp, f_side1, args.win)
    f_chr2, f_s2, f_e2 = bp_to_region(f_chr2_bp, f_pos2_bp, f_side2, args.win)

    # 统计真结构 / 假结构 read-pairs
    true_pairs = count_pairs_between_regions(
        bam, t_chr1, t_s1, t_e1, t_chr2, t_s2, t_e2
    )
    false_pairs = count_pairs_between_regions(
        bam, f_chr1, f_s1, f_e1, f_chr2, f_s2, f_e2
    )

    denom = true_pairs + false_pairs
    if denom > 0:
        true_ratio = true_pairs / denom
        false_ratio = false_pairs / denom
    else:
        true_ratio = float("nan")
        false_ratio = float("nan")

    header = "\t".join([
        "bam",
        "win_bp",
        "true_bp1", "true_bp2",
        "true_region1", "true_region2", "true_pairs",
        "false_bp1", "false_bp2",
        "false_region1", "false_region2", "false_pairs",
        "true_ratio", "false_ratio"
    ]) + "\n"

    true_region1 = f"{t_chr1}:{t_s1}-{t_e1}"
    true_region2 = f"{t_chr2}:{t_s2}-{t_e2}"
    false_region1 = f"{f_chr1}:{f_s1}-{f_e1}"
    false_region2 = f"{f_chr2}:{f_s2}-{f_e2}"

    line = "\t".join([
        args.bam,
        str(args.win),
        args.true_bp1, args.true_bp2,
        true_region1, true_region2, str(true_pairs),
        args.false_bp1, args.false_bp2,
        false_region1, false_region2, str(false_pairs),
        f"{true_ratio:.4f}", f"{false_ratio:.4f}",
    ]) + "\n"

    if args.out == "-" or args.out == "stdout":
        if not args.no_header:
            print(header, end="")
        print(line, end="")
    else:
        with open(args.out, "a") as f:
            if not args.no_header:
                f.write(header)
            f.write(line)

    bam.close()

if __name__ == "__main__":
    main()
