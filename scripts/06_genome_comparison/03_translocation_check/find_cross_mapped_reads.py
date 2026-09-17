#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pysam
from collections import Counter


# ----------------- 基础工具函数 -----------------

def parse_region(s):
    """
    解析形如 Chr1:1-22839936 的字符串
    返回 (chrom, start, end) （1-based, 闭区间）
    """
    chrom, coords = s.split(":")
    start, end = coords.split("-")
    return chrom, int(start), int(end)


def in_region(chrom, pos1based, target_chr, start, end):
    """判断 (chrom, pos) 是否落在指定区域内"""
    return (chrom == target_chr) and (start <= pos1based <= end)


def collect_read1_names_in_region(bam_path, region_str):
    """
    在给定 BAM + 区域里收集所有落在区域内的 read1 名称集合
    只看 primary、成对、两端都 mapped
    """
    chrom, start, end = parse_region(region_str)
    names = set()

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        start0 = max(start - 1, 0)
        end0 = end
        for read in bam.fetch(chrom, start0, end0):
            if not read.is_paired or not read.is_read1:
                continue
            if read.is_unmapped or read.mate_is_unmapped:
                continue
            if read.is_secondary or read.is_supplementary:
                continue

            pos1based = read.reference_start + 1
            if in_region(chrom, pos1based, chrom, start, end):
                names.add(read.query_name)

    return names


def summarize_mate_positions(bam_path, region_str, read_names):
    """
    对于给定 BAM 和 区域 + read_name 集合：
    - 找到所有 read1 落在该区域且 read_name 属于 read_names
    - 统计 mate 比对到的染色体数量
    - 返回:
        mate_chr_counter: Counter(chr -> count)
        mate_records: list[(qname, read_chr, read_pos, mate_chr, mate_pos)]
    """
    chrom, start, end = parse_region(region_str)
    mate_chr_counter = Counter()
    mate_records = []

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        start0 = max(start - 1, 0)
        end0 = end
        for read in bam.fetch(chrom, start0, end0):
            if not read.is_paired or not read.is_read1:
                continue
            if read.is_unmapped or read.mate_is_unmapped:
                continue
            if read.is_secondary or read.is_supplementary:
                continue

            if read.query_name not in read_names:
                continue

            rchr = bam.get_reference_name(read.reference_id)
            rpos = read.reference_start + 1

            mate_id = read.next_reference_id
            if mate_id < 0:
                continue
            mchr = bam.get_reference_name(mate_id)
            mpos = read.next_reference_start + 1

            mate_chr_counter[mchr] += 1
            mate_records.append((read.query_name, rchr, rpos, mchr, mpos))

    return mate_chr_counter, mate_records


def write_subset_bam(in_bam_path, out_bam_path, read_name_set):
    """
    从原始 BAM 中抽取 query_name 属于 read_name_set 的所有 reads（两端都写出）
    """
    with pysam.AlignmentFile(in_bam_path, "rb") as bam_in:
        with pysam.AlignmentFile(out_bam_path, "wb", template=bam_in) as bam_out:
            for read in bam_in.fetch(until_eof=True):
                if read.query_name in read_name_set:
                    bam_out.write(read)


def fetch_positions_from_bam(bam_path, read_names):
    """
    在 BAM 中查找这些 read_names 的 read1 比对位置 + mate 位置
    返回字典: name -> (read_chr, read_pos, mate_chr, mate_pos)
    """
    name_set = set(read_names)
    result = {}

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.query_name not in name_set:
                continue
            if not read.is_paired or not read.is_read1:
                continue
            if read.is_unmapped or read.mate_is_unmapped:
                continue
            if read.is_secondary or read.is_supplementary:
                continue

            rchr = bam.get_reference_name(read.reference_id)
            rpos = read.reference_start + 1
            mate_id = read.next_reference_id
            if mate_id < 0:
                continue
            mchr = bam.get_reference_name(mate_id)
            mpos = read.next_reference_start + 1

            result[read.query_name] = (rchr, rpos, mchr, mpos)

            if len(result) == len(name_set):
                break

    return result


def pair_distance(chr_a, pos_a, chr_b, pos_b):
    """
    计算一对 reads 之间的距离：
    - 如果在同一条染色体: abs(pos_b - pos_a)
    - 如果跨染色体: 返回 -1 （可当 NA）
    """
    if chr_a == chr_b:
        return abs(pos_b - pos_a)
    else:
        return -1


# ----------------- 主程序 -----------------

def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "同一套 reads 比对到两个基因组："
            "在 bam-gifu 中落在 region-gifu，且在 bam-mg 中落在 region-mg 的 read1 交集；\n"
            "在两边分别找出 Chr1<->Chr2 的跨染色体 pairs，"
            "输出 merged 表（两基因组上的位置 + pair 距离），"
            "并导出对应的小 BAM。"
        )
    )
    p.add_argument("--bam-gifu", required=True,
                   help="同一套 reads 比对到 Gifu 基因组的 BAM")
    p.add_argument("--bam-mg", required=True,
                   help="同一套 reads 比对到 MG20 基因组的 BAM")
    p.add_argument("--region-gifu", required=True,
                   help="Gifu 区域，如 Chr1:22829936-22839936")
    p.add_argument("--region-mg", required=True,
                   help="MG20 区域，如 Chr2:99270160-99280160")
    p.add_argument("--out-prefix", default="cross_reads",
                   help="输出文件前缀")
    return p.parse_args()


def main():
    args = parse_args()

    # 把 region 中的染色体名当作“Chr1”和“Chr2”的标签
    gifu_chr, _, _ = parse_region(args.region_gifu)   # 例如 Chr1
    mg_chr, _, _   = parse_region(args.region_mg)     # 例如 Chr2

    chr1_name = gifu_chr
    chr2_name = mg_chr

    # ---------- 1. 区域内 read1 名称 & 交集 ----------
    print("[*] 收集 Gifu BAM 区域内的 read1 名称...")
    names_gifu = collect_read1_names_in_region(args.bam_gifu, args.region_gifu)
    print(f"    Gifu 区域内 read1 数量: {len(names_gifu)}")

    print("[*] 收集 MG20 BAM 区域内的 read1 名称...")
    names_mg = collect_read1_names_in_region(args.bam_mg, args.region_mg)
    print(f"    MG20 区域内 read1 数量: {len(names_mg)}")

    cross_names = names_gifu & names_mg
    print(f"[*] 同时满足两个区域的 read1 数量（交集）: {len(cross_names)}")
    if not cross_names:
        print("[!] 交集为空，结束。")
        return

    # ---------- 2. 在两边 BAM 中统计 mate 位置 ----------
    print("[*] 分析这些 reads 在 Gifu BAM 中的 mate 位置...")
    gifu_mate_chr_counter, gifu_mate_records = summarize_mate_positions(
        args.bam_gifu, args.region_gifu, cross_names
    )

    print("[*] 分析这些 reads 在 MG20 BAM 中的 mate 位置...")
    mg_mate_chr_counter, mg_mate_records = summarize_mate_positions(
        args.bam_mg, args.region_mg, cross_names
    )

    print("\n=== Summary: mate 染色体分布（Gifu BAM） ===")
    for chr_name, count in gifu_mate_chr_counter.most_common():
        print(f"Gifu_BAM_mate_chr\t{chr_name}\t{count}")

    print("\n=== Summary: mate 染色体分布（MG20 BAM） ===")
    for chr_name, count in mg_mate_chr_counter.most_common():
        print(f"MG20_BAM_mate_chr\t{chr_name}\t{count}")

    # ---------- 3. 各自 BAM 中的 Chr1<->Chr2 跨染色体 pairs ----------
    gifu_cross_records = []
    gifu_cross_names = set()
    for qname, rchr, rpos, mchr, mpos in gifu_mate_records:
        if ((rchr == chr1_name and mchr == chr2_name) or
            (rchr == chr2_name and mchr == chr1_name)):
            gifu_cross_records.append((qname, rchr, rpos, mchr, mpos))
            gifu_cross_names.add(qname)

    mg_cross_records = []
    mg_cross_names = set()
    for qname, rchr, rpos, mchr, mpos in mg_mate_records:
        if ((rchr == chr1_name and mchr == chr2_name) or
            (rchr == chr2_name and mchr == chr1_name)):
            mg_cross_records.append((qname, rchr, rpos, mchr, mpos))
            mg_cross_names.add(qname)

    print("\n=== Chr1<->Chr2 跨染色体 pairs 数量 ===")
    print(f"Gifu BAM 中: {len(gifu_cross_records)} 条")
    print(f"MG20 BAM 中: {len(mg_cross_records)} 条")

    # ---------- 4A. 以 Gifu 为锚点的 merged 表 ----------
    if gifu_cross_records:
        print("[*] 生成 Gifu-anchored merged 表...")
        mg_pos_map_for_gifu = fetch_positions_from_bam(args.bam_mg, gifu_cross_names)

        merged_tsv_gifu = f"{args.out_prefix}.gifu_anchored.merged_pairs.tsv"
        with open(merged_tsv_gifu, "w") as out:
            out.write(
                "read_name\t"
                "gifu_read_chr\tgifu_read_pos\tgifu_mate_chr\tgifu_mate_pos\tgifu_pair_dist\t"
                "mg_read_chr\tmg_read_pos\tmg_mate_chr\tmg_mate_pos\tmg_pair_dist\n"
            )
            for qname, g_rchr, g_rpos, g_mchr, g_mpos in gifu_cross_records:
                g_dist = pair_distance(g_rchr, g_rpos, g_mchr, g_mpos)
                if qname in mg_pos_map_for_gifu:
                    m_rchr, m_rpos, m_mchr, m_mpos = mg_pos_map_for_gifu[qname]
                    m_dist = pair_distance(m_rchr, m_rpos, m_mchr, m_mpos)
                else:
                    m_rchr = m_mchr = "NA"
                    m_rpos = m_mpos = -1
                    m_dist = -1
                out.write(
                    f"{qname}\t"
                    f"{g_rchr}\t{g_rpos}\t{g_mchr}\t{g_mpos}\t{g_dist}\t"
                    f"{m_rchr}\t{m_rpos}\t{m_mchr}\t{m_mpos}\t{m_dist}\n"
                )
        print(f"    写入: {merged_tsv_gifu}")

        # 导出对应的小 BAM
        gifu_pairs_bam = f"{args.out_prefix}.gifu_chr1_chr2_pairs.bam"
        print(f"[*] 从 {args.bam_gifu} 中抽取 {len(gifu_cross_names)} 条 reads 写入 {gifu_pairs_bam}")
        write_subset_bam(args.bam_gifu, gifu_pairs_bam, gifu_cross_names)

        gifu_reads_in_mg_bam = f"{args.out_prefix}.gifu_reads_in_mg.bam"
        print(f"[*] 从 {args.bam_mg} 中抽取同名 reads 写入 {gifu_reads_in_mg_bam}")
        write_subset_bam(args.bam_mg, gifu_reads_in_mg_bam, gifu_cross_names)

    # ---------- 4B. 以 MG 为锚点的 merged 表（一般你更关注这一部分） ----------
    if mg_cross_records:
        print("[*] 生成 MG-anchored merged 表...")
        gifu_pos_map_for_mg = fetch_positions_from_bam(args.bam_gifu, mg_cross_names)

        merged_tsv_mg = f"{args.out_prefix}.mg_anchored.merged_pairs.tsv"
        with open(merged_tsv_mg, "w") as out:
            out.write(
                "read_name\t"
                "mg_read_chr\tmg_read_pos\tmg_mate_chr\tmg_mate_pos\tmg_pair_dist\t"
                "gifu_read_chr\tgifu_read_pos\tgifu_mate_chr\tgifu_mate_pos\tgifu_pair_dist\n"
            )
            for qname, m_rchr, m_rpos, m_mchr, m_mpos in mg_cross_records:
                m_dist = pair_distance(m_rchr, m_rpos, m_mchr, m_mpos)
                if qname in gifu_pos_map_for_mg:
                    g_rchr, g_rpos, g_mchr, g_mpos = gifu_pos_map_for_mg[qname]
                    g_dist = pair_distance(g_rchr, g_rpos, g_mchr, g_mpos)
                else:
                    g_rchr = g_mchr = "NA"
                    g_rpos = g_mpos = -1
                    g_dist = -1
                out.write(
                    f"{qname}\t"
                    f"{m_rchr}\t{m_rpos}\t{m_mchr}\t{m_mpos}\t{m_dist}\t"
                    f"{g_rchr}\t{g_rpos}\t{g_mchr}\t{g_mpos}\t{g_dist}\n"
                )
        print(f"    写入: {merged_tsv_mg}")

        # 导出对应的小 BAM
        mg_pairs_bam = f"{args.out_prefix}.mg_chr1_chr2_pairs.bam"
        print(f"[*] 从 {args.bam_mg} 中抽取 {len(mg_cross_names)} 条 reads 写入 {mg_pairs_bam}")
        write_subset_bam(args.bam_mg, mg_pairs_bam, mg_cross_names)

        mg_reads_in_gifu_bam = f"{args.out_prefix}.mg_reads_in_gifu.bam"
        print(f"[*] 从 {args.bam_gifu} 中抽取同名 reads 写入 {mg_reads_in_gifu_bam}")
        write_subset_bam(args.bam_gifu, mg_reads_in_gifu_bam, mg_cross_names)


if __name__ == "__main__":
    main()