#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
批量调用 find_cross_mapped_reads.py，对多个材料做“Gifu 型 / MG20 型”结构分型。

输入：
  一个 3 列的配置文件，比如 trans.conf.txt：
      sample   bam_gifu   bam_mg
  或你当前用的 lotus_translocation_typing.conf 那种格式。

含义：
  - bam_gifu : 该样本 reads 比对到 Gifu 基因组的 BAM
  - bam_mg   : 该样本 reads 比对到 MG20 基因组的 BAM

对于每个样本：
  调用 find_cross_mapped_reads.py 生成：
    {sample}.trans.gifu_anchored.merged_pairs.tsv
    {sample}.trans.mg_anchored.merged_pairs.tsv

  统计：
    n_gifu = gifu_anchored.merged_pairs.tsv 行数（在 Gifu 参考上 Chr1↔Chr2 跨染色体 pairs）
    n_mg   = mg_anchored.merged_pairs.tsv 行数（在 MG20 参考上 Chr1↔Chr2 跨染色体 pairs）

  判定逻辑（按你的要求）：
    - 如果 n_mg ≥ min_support 且 n_gifu < min_support  → Gifu-like
    - 如果 n_gifu ≥ min_support 且 n_mg < min_support  → MG20-like
    - 如果 n_gifu ≥ min_support 且 n_mg ≥ min_support → Mixed
    - 否则                                         → Uncertain

  默认 min_support = 5 （≥5 条就可以判断）
"""

import argparse
import subprocess
import os
import pandas as pd


def run_find_cross(bam_gifu, bam_mg, region_gifu, region_mg,
                   out_prefix, script_path):
    """调用 find_cross_mapped_reads.py"""
    cmd = [
        "python", script_path,
        "--bam-gifu", bam_gifu,
        "--bam-mg", bam_mg,
        "--region-gifu", region_gifu,
        "--region-mg", region_mg,
        "--out-prefix", out_prefix,
    ]
    print("[CMD]", " ".join(cmd))
    subprocess.run(cmd, check=True)


def count_rows_if_exists(path):
    """存在就读 TSV 数行，不存在返回 0"""
    if os.path.exists(path):
        df = pd.read_csv(path, sep="\t")
        return df.shape[0]
    else:
        return 0


def classify_sample(n_gifu, n_mg, min_support=5):
    """
    根据 GifuBAM / MG20BAM 中跨染色体 pair 数量给出分型（方向已经按你的要求反过来）：

      - n_mg >= min_support 且 n_gifu < min_support  → Gifu-like
      - n_gifu >= min_support 且 n_mg < min_support → MG20-like
      - 两边都 >= min_support                       → Mixed
      - 否则                                         → Uncertain
    """
    total = n_gifu + n_mg
    ratio_mg = (n_mg / total) if total > 0 else 0.0

    # 只在 MG20 参考中有跨染色体 pair：说明样本本身结构更像 Gifu
    if n_mg >= min_support and n_gifu < min_support:
        return "Gifu-like", ratio_mg, "cross only in MG20-BAM"

    # 只在 Gifu 参考中有跨染色体 pair：说明样本本身结构更像 MG20
    if n_gifu >= min_support and n_mg < min_support:
        return "MG20-like", ratio_mg, "cross only in Gifu-BAM"

    # 两边都有足够的跨染色体 pair：可能是混合群体 / 杂合结构
    if n_gifu >= min_support and n_mg >= min_support:
        return "Mixed", ratio_mg, "cross in both BAMs"

    # 支持数太少：不确定
    return "Uncertain", ratio_mg, f"support< {min_support} on both sides"


def main():
    ap = argparse.ArgumentParser(
        description="Batch typing of Chr1↔Chr2 translocation (Gifu vs MG20)"
    )
    ap.add_argument("--conf", required=True,
                    help="3 列: sample  bam_gifu  bam_mg")
    ap.add_argument("--region-gifu", required=True,
                    help="Gifu 断点区域，如 Chr1:22829936-22839936")
    ap.add_argument("--region-mg", required=True,
                    help="MG20 断点区域，如 Chr2:99270160-99280160")
    ap.add_argument("--find-script", default="find_cross_mapped_reads.py",
                    help="find_cross_mapped_reads.py 的路径")
    ap.add_argument("--out-summary", default="translocation_typing.summary.tsv",
                    help="输出 summary 表")
    ap.add_argument("--min-support", type=int, default=5,
                    help="最少跨染色体 pairs 支持数（默认 5）")
    args = ap.parse_args()

    # 读取配置文件
    conf = pd.read_csv(args.conf, sep="\t", header=None,
                       names=["sample", "bam_gifu", "bam_mg"])

    records = []

    for _, row in conf.iterrows():
        sample = row["sample"]
        bam_gifu = row["bam_gifu"]
        bam_mg = row["bam_mg"]

        out_prefix = f"{sample}.trans"
        print(f"\n===== Processing sample: {sample} =====")
        print(f"  bam_gifu: {bam_gifu}")
        print(f"  bam_mg  : {bam_mg}")

        # 1. 调用 find_cross_mapped_reads.py
        run_find_cross(
            bam_gifu=bam_gifu,
            bam_mg=bam_mg,
            region_gifu=args.region_gifu,
            region_mg=args.region_mg,
            out_prefix=out_prefix,
            script_path=args.find_script,
        )

        # 2. 统计 gifu-anchored / mg-anchored cross pairs
        gifu_tsv = f"{out_prefix}.gifu_anchored.merged_pairs.tsv"
        mg_tsv = f"{out_prefix}.mg_anchored.merged_pairs.tsv"

        n_gifu = count_rows_if_exists(gifu_tsv)
        n_mg = count_rows_if_exists(mg_tsv)

        call, ratio_mg, note = classify_sample(
            n_gifu=n_gifu,
            n_mg=n_mg,
            min_support=args.min_support,
        )

        print(f"  n_cross_in_gifu_bam = {n_gifu}")
        print(f"  n_cross_in_mg_bam   = {n_mg}")
        print(f"  => call = {call}  (ratio_MG20 = {ratio_mg:.3f})  note: {note}")

        records.append({
            "sample": sample,
            "n_cross_in_gifu_bam": n_gifu,
            "n_cross_in_mg_bam": n_mg,
            "call": call,
            "ratio_MG20": f"{ratio_mg:.4f}",
            "note": note,
        })

    summary_df = pd.DataFrame(records)
    summary_df.to_csv(args.out_summary, sep="\t", index=False)
    print(f"\n[Done] Summary written to: {args.out_summary}")


if __name__ == "__main__":
    main()
