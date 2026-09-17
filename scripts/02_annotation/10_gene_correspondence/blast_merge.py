#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import pandas as pd
from collections import defaultdict


def load_length_map(length_file: str, key_col: str, val_col: str) -> dict:
    df = pd.read_csv(length_file, sep="\t", header=None, usecols=[0, 1])
    df.columns = [key_col, val_col]
    df[key_col] = df[key_col].astype(str)
    df[val_col] = df[val_col].astype(float)
    return dict(zip(df[key_col], df[val_col]))


def process_blast_results(input_file, query_lengths_file, reference_lengths_file, output_file,
                          chunksize=2_000_000):

    qlen_map = load_length_map(query_lengths_file, "query_id", "query_length")
    rlen_map = load_length_map(reference_lengths_file, "reference_id", "reference_length")

    # (q, r) -> 累加值
    acc_sim_len = defaultdict(float)   # Σ(similarity * aln_len)
    acc_len = defaultdict(float)       # Σ(aln_len)

    dtypes = {0: "string", 1: "string", 2: "float32", 3: "float32"}

    reader = pd.read_csv(
        input_file,
        sep="\t",
        header=None,
        usecols=[0, 1, 2, 3],
        dtype=dtypes,
        chunksize=chunksize
    )

    for chunk in reader:
        chunk.columns = ["query_id", "reference_id", "similarity", "alignment_length"]
        chunk["sim_x_len"] = chunk["similarity"] * chunk["alignment_length"]

        g = chunk.groupby(["query_id", "reference_id"], sort=False, observed=True).agg(
            sum_sim_x_len=("sim_x_len", "sum"),
            sum_len=("alignment_length", "sum")
        )

        for (q, r), row in g.iterrows():
            key = (str(q), str(r))
            acc_sim_len[key] += float(row["sum_sim_x_len"])
            acc_len[key] += float(row["sum_len"])

    # ===== 输出 =====
    with open(output_file, "w") as out:
        for (q, r), total_len in acc_len.items():
            if total_len <= 0:
                continue

            sim_sum = acc_sim_len[(q, r)]
            sim = sim_sum / total_len

            qlen = qlen_map.get(q, float("nan"))
            rlen = rlen_map.get(r, float("nan"))

            qcov = total_len / qlen if qlen == qlen else float("nan")
            scov = total_len / rlen if rlen == rlen else float("nan")

            # ===== 新增两列 =====
            q_weighted_sim_cov = sim_sum / qlen if qlen == qlen else float("nan")
            r_weighted_sim_cov = sim_sum / rlen if rlen == rlen else float("nan")

            out.write(
                f"{q}\t{r}\t{sim:.6f}\t{total_len:.2f}\t"
                f"{qlen if qlen == qlen else 'NA'}\t"
                f"{rlen if rlen == rlen else 'NA'}\t"
                f"{qcov if qcov == qcov else 'NA'}\t"
                f"{scov if scov == scov else 'NA'}\t"
                f"{q_weighted_sim_cov if q_weighted_sim_cov == q_weighted_sim_cov else 'NA'}\t"
                f"{r_weighted_sim_cov if r_weighted_sim_cov == r_weighted_sim_cov else 'NA'}\n"
            )


if __name__ == "__main__":
    if len(sys.argv) not in (5, 6):
        print(
            "用法: python3 blast_merge.py "
            "<blast_result> <query_length_file> <ref_length_file> <output_file> [chunksize]"
        )
        sys.exit(1)

    input_file = sys.argv[1]
    query_lengths_file = sys.argv[2]
    reference_lengths_file = sys.argv[3]
    output_file = sys.argv[4]
    chunksize = int(sys.argv[5]) if len(sys.argv) == 6 else 2_000_000

    process_blast_results(
        input_file,
        query_lengths_file,
        reference_lengths_file,
        output_file,
        chunksize=chunksize
    )
