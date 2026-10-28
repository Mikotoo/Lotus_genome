#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
from collections import defaultdict, Counter

LEN_RE = re.compile(r"\|L=(\d+)\|")

def parse_len(seq_id: str) -> int:
    m = LEN_RE.search(seq_id)
    if not m:
        raise ValueError(f"Cannot parse length from id: {seq_id}")
    return int(m.group(1))

def score_candidate(L, lengths, tol=6, max_k=6):
    """
    Score how well a candidate monomer length L explains the observed lengths.
    A length x is considered explained if exists k s.t. |x - k*L| <= tol.
    Return: (explained_count, explained_weighted, details)
    """
    explained = 0
    explained_w = 0.0
    detail = []
    for x in lengths:
        best = None
        for k in range(1, max_k + 1):
            target = k * L
            d = abs(x - target)
            if d <= tol:
                if best is None or d < best[0]:
                    best = (d, k, target)
        if best:
            explained += 1
            # closer match contributes more
            explained_w += 1.0 / (1.0 + best[0])
            detail.append((x, best[1], best[2], best[0]))
        else:
            detail.append((x, None, None, None))
    return explained, explained_w, detail

def main():
    ap = argparse.ArgumentParser(
        description="Pick shortest monomer per cluster by explaining k×L length structure."
    )
    ap.add_argument("-i", "--input", required=True, help="clusters.tsv (cluster_id\\tmembers...)")
    ap.add_argument("-o", "--output", default="cluster_monomer.tsv", help="Output TSV")
    ap.add_argument("--tol", type=int, default=6, help="Tolerance in bp for k×L matching (default: 6)")
    ap.add_argument("--max-k", type=int, default=6, help="Max multiplier k to try (default: 6)")
    ap.add_argument("--top", type=int, default=5, help="Report top N monomer-length member IDs closest to L* (default: 5)")
    args = ap.parse_args()

    # cluster -> list of member IDs
    clusters = {}
    with open(args.input, "r", encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        # accept either: cluster_id size members  OR your pasted format cluster_id members...
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            cid = parts[0]
            if len(parts) == 2:
                members = parts[1].split(",") if "," in parts[1] else parts[1].split()
            elif len(parts) >= 3 and parts[1].isdigit():
                members = parts[2].split(",")
            else:
                # fallback: everything after first col is member IDs (space separated)
                members = []
                for p in parts[1:]:
                    members.extend(p.split(","))
            members = [m for m in members if m]
            clusters[cid] = members

    with open(args.output, "w", encoding="utf-8") as out:
        out.write("cluster\tbest_monomer_len\tn_members\tn_distinct_lens\texplained_count\texplained_weighted\tlength_histogram\ttop_monomer_like_ids\n")

        for cid, members in sorted(clusters.items()):
            lens = [parse_len(m) for m in members]
            lens_sorted = sorted(lens)
            hist = Counter(lens_sorted)

            # candidate monomer lengths: use observed lengths in a reasonable range
            # for centromere satellites, monomer often 80–250, but keep broad
            candidates = sorted({L for L in hist if 60 <= L <= 250})

            if not candidates:
                # fallback: take min observed as candidate
                candidates = [min(lens_sorted)]

            best = None
            best_detail = None
            for L in candidates:
                expl, expl_w, detail = score_candidate(L, lens_sorted, tol=args.tol, max_k=args.max_k)
                # primary: maximize explained_count, secondary: explained_weighted, tertiary: smaller L
                key = (expl, expl_w, -L)
                if best is None or key > best[0]:
                    best = (key, L, expl, expl_w)
                    best_detail = detail

            _, bestL, expl, expl_w = best

            # pick member IDs closest to bestL (k=1)
            # choose those with length within tol of bestL, or just closest
            monomer_like = [(abs(parse_len(m) - bestL), parse_len(m), m) for m in members]
            monomer_like.sort(key=lambda x: (x[0], x[1]))
            top_ids = [x[2] for x in monomer_like[:args.top]]

            # histogram string
            hist_str = ",".join([f"{L}:{c}" for L, c in sorted(hist.items())])

            out.write(f"{cid}\t{bestL}\t{len(members)}\t{len(hist)}\t{expl}\t{expl_w:.3f}\t{hist_str}\t" + ",".join(top_ids) + "\n")

    print(f"[DONE] Wrote: {args.output}")

if __name__ == "__main__":
    main()
