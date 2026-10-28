#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
from pathlib import Path
import argparse
import sys

EPS = 1e-6

def calc_tau(row: pd.Series) -> float:
    """
    τ (tau) tissue-specificity index.
    row: expression across tissues (non-negative)
    """
    xmax = row.max()
    if xmax <= 0:
        return 0.0
    n = row.shape[0]
    return float(((1.0 - row / xmax).sum()) / (n - 1))

def build_tissue_matrix(df: pd.DataFrame, tissue_patterns: dict) -> pd.DataFrame:
    """
    df: gene_id as index, columns are sample names (TPM)
    tissue_patterns: {tissue: [substr1, substr2, ...]}  (OR match)
    returns: DataFrame (gene_id x tissue) with mean TPM across matched columns
    """
    tissue_expr = {}
    for tissue, pats in tissue_patterns.items():
        cols = [c for c in df.columns if any(p in c for p in pats)]
        if len(cols) == 0:
            raise ValueError(f"[ERROR] No columns matched for tissue '{tissue}' using patterns: {pats}")
        tissue_expr[tissue] = df[cols].mean(axis=1)
    return pd.DataFrame(tissue_expr)

def process_one(expr_tsv: str, label: str, outdir: str, topn: int,
                tissue_patterns: dict, min_tpm_max: float = 0.0) -> dict:
    """
    Read TPM matrix, compute tissue mean, tau, max_tissue, and export topN lists per tissue.
    Returns: dict {tissue: set(gene_id)}
    """
    df = pd.read_csv(expr_tsv, sep="\t", dtype={0: str})
    if "gene_id" not in df.columns:
        raise ValueError(f"[ERROR] '{expr_tsv}' has no column named 'gene_id'")
    df = df.set_index("gene_id")

    # Ensure numeric
    df = df.apply(pd.to_numeric, errors="coerce").fillna(0.0)

    # 1) tissue-level mean TPM
    tissue_df = build_tissue_matrix(df, tissue_patterns)

    # Optional filter: only keep genes expressed above some level in at least one tissue
    if min_tpm_max > 0:
        keep = tissue_df.max(axis=1) >= min_tpm_max
        tissue_df = tissue_df.loc[keep]

    # 2) tau + max tissue
    tau = tissue_df.apply(calc_tau, axis=1)
    max_tissue = tissue_df.idxmax(axis=1)

    res = tissue_df.copy()
    res["tau"] = tau
    res["max_tissue"] = max_tissue

    # 3) export topN per tissue
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    top_sets = {}
    for tissue in tissue_patterns.keys():
        sub = res[res["max_tissue"] == tissue].sort_values("tau", ascending=False)
        sub_top = sub.head(topn)
        genes = list(sub_top.index)
        top_sets[tissue] = set(genes)

        # write list
        (out / f"{label}_{tissue}_top{topn}.txt").write_text("\n".join(genes) + "\n", encoding="utf-8")

        # optional: write a table for debugging/traceability
        sub_top.to_csv(out / f"{label}_{tissue}_top{topn}.detail.tsv", sep="\t", index=True)

    # write a full summary table (tissue means + tau + max_tissue)
    res.to_csv(out / f"{label}.tissue_mean_tau.tsv", sep="\t", index=True)

    return top_sets

def jaccard(a: set, b: set) -> float:
    if len(a) == 0 and len(b) == 0:
        return 1.0
    if len(a) == 0 or len(b) == 0:
        return 0.0
    return len(a & b) / len(a | b)

def main():
    ap = argparse.ArgumentParser(
        description="Compute tau tissue-specificity, pick topN tissue-enriched genes, and compute Jaccard between Gifu and MG20."
    )
    ap.add_argument("--gifu", required=True, help="Gifu TPM matrix TSV (gene_id + samples)")
    ap.add_argument("--mg20", required=True, help="MG20 TPM matrix TSV (gene_id + samples)")
    ap.add_argument("--topn", type=int, default=500, help="Top N genes per tissue (default: 500)")
    ap.add_argument("--min_tpm_max", type=float, default=0.0,
                    help="Optional filter: keep genes with max tissue-mean TPM >= this value (default: 0, no filter)")
    ap.add_argument("--outprefix", default="tau_topN_jaccard", help="Output prefix (default: tau_topN_jaccard)")
    args = ap.parse_args()

    # ========= Tissue definitions (root_uni is separated) =========
    # NOTE: patterns are substrings; OR logic within each tissue
    tissue_patterns = {
        "flower":   ["flower_"],
        "leaf":     ["leaf_"],
        "nodule":   ["nodule_"],
        "pod":      ["pod_"],
        "stem":     ["stem_"],

        # root excludes root_uni explicitly by using dpi/hpi patterns
        "root":     ["root_dpi", "root_hpi"],

        # root_uni is its own group
        "root_uni": ["root_uni_"],
    }

    gifu_out = f"{args.outprefix}.Gifu"
    mg20_out = f"{args.outprefix}.MG20"

    print(f"[INFO] Processing Gifu: {args.gifu}")
    gifu_sets = process_one(
        args.gifu, "Gifu", gifu_out, args.topn, tissue_patterns, min_tpm_max=args.min_tpm_max
    )

    print(f"[INFO] Processing MG20: {args.mg20}")
    mg20_sets = process_one(
        args.mg20, "MG20", mg20_out, args.topn, tissue_patterns, min_tpm_max=args.min_tpm_max
    )

    # ========= Jaccard per tissue =========
    rows = []
    for tissue in tissue_patterns.keys():
        A = gifu_sets.get(tissue, set())
        B = mg20_sets.get(tissue, set())
        inter = len(A & B)
        union = len(A | B)
        jac = jaccard(A, B)
        rows.append([tissue, len(A), len(B), inter, union, jac])

    out = pd.DataFrame(rows, columns=["tissue", "Gifu_topN", "MG20_topN", "intersection", "union", "jaccard"])
    out_tsv = f"{args.outprefix}.jaccard_by_tissue.tsv"
    out.to_csv(out_tsv, sep="\t", index=False)

    print("[OK] Jaccard results:")
    print(out.to_string(index=False))
    print(f"[OUT] {out_tsv}")
    print(f"[OUT] TopN lists and details: {gifu_out}/ , {mg20_out}/")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(str(e), file=sys.stderr)
        sys.exit(1)
