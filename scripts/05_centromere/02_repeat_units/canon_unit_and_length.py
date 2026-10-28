#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from collections import defaultdict

COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")

def revcomp(seq: str) -> str:
    return seq.translate(COMP)[::-1]

def canonical_rotation(seq: str) -> str:
    """Return lexicographically smallest rotation of seq."""
    if not seq:
        return seq
    # O(n^2) but n is small (<~200) for satellite monomers, ok.
    rots = (seq[i:] + seq[:i] for i in range(len(seq)))
    return min(rots)

def canonical_unit_strand_invariant(seq: str) -> str:
    """Canonicalize by rotation, and also consider reverse-complement."""
    s = seq.upper()
    c1 = canonical_rotation(s)
    rc = revcomp(s)
    c2 = canonical_rotation(rc)
    return min(c1, c2)

def is_header(line: str) -> bool:
    # ignore TRF header-like lines if present
    low = line.lower()
    return (low.startswith("sequence:") or low.startswith("parameters:"))

def main():
    ap = argparse.ArgumentParser(
        description="Compute canonical unit (rotation + reverse-complement) and unit length from TRF-derived table."
    )
    ap.add_argument("-i", "--input", required=True, help="Input TSV (or whitespace-separated) file.")
    ap.add_argument("-o", "--output", required=True, help="Output TSV with canonical unit columns added.")
    ap.add_argument("--seq-col", type=int, default=4,
                    help="1-based column index for unit sequence (default: 4).")
    ap.add_argument("--start-col", type=int, default=2,
                    help="1-based column index for start (default: 2).")
    ap.add_argument("--end-col", type=int, default=3,
                    help="1-based column index for end (default: 3).")
    ap.add_argument("--summary", default="length_summary.tsv",
                    help="Output summary TSV filename (default: length_summary.tsv).")
    args = ap.parse_args()

    seq_idx = args.seq_col - 1
    s_idx = args.start_col - 1
    e_idx = args.end_col - 1

    # length -> count, total_span
    len_count = defaultdict(int)
    len_span = defaultdict(int)

    with open(args.input, "r", encoding="utf-8") as fin, \
         open(args.output, "w", encoding="utf-8") as fout:

        # 写表头
        fout.write("\t".join([
            "chrom", "start", "end", "unit",
            "unit_len", "canon_unit", "canon_len",
            "span_bp"
        ]) + "\t" + "raw_line" + "\n")

        for line in fin:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            if is_header(line):
                continue

            # 允许 TSV 或空格分隔
            parts = line.split("\t")
            if len(parts) == 1:
                parts = line.split()

            # 跳过明显不是数据行的情况
            if len(parts) <= max(seq_idx, s_idx, e_idx):
                continue

            chrom = parts[0]
            try:
                start = int(parts[s_idx])
                end = int(parts[e_idx])
            except ValueError:
                # 可能是表头行
                continue

            unit = parts[seq_idx].strip()
            # 若 unit 中混入了非ACGTN字符（一般不会），这里简单清理
            unit_clean = "".join([c for c in unit.upper() if c in "ACGTN"])
            if not unit_clean:
                continue

            unit_len = len(unit_clean)
            canon = canonical_unit_strand_invariant(unit_clean)
            span = end - start + 1
            if span < 0:
                span = 0

            len_count[unit_len] += 1
            len_span[unit_len] += span

            # 输出：固定字段 + 原始整行（方便你回溯）
            fout.write("\t".join(map(str, [
                chrom, start, end, unit_clean,
                unit_len, canon, len(canon),
                span
            ])) + "\t" + line + "\n")

    # 输出长度统计
    with open(args.summary, "w", encoding="utf-8") as fs:
        fs.write("unit_len\tcount\ttotal_span_bp\n")
        for L in sorted(len_count.keys()):
            fs.write(f"{L}\t{len_count[L]}\t{len_span[L]}\n")

if __name__ == "__main__":
    main()
