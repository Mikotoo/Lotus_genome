#!/usr/bin/env python3
import argparse
from pathlib import Path


def parse_args():
    ap = argparse.ArgumentParser(
        description="Build superclusters with Union-Find from a pairs table (A,B,Jaccard)."
    )
    ap.add_argument("--outdir", required=True, help="Working directory")
    ap.add_argument(
        "--pairs",
        default=None,
        help="Pairs TSV path. Default: <outdir>/pairs.jaccard_ge<JAC>.tsv if --jac is given, else <outdir>/pairs.jaccard.tsv",
    )
    ap.add_argument(
        "--clusters",
        default=None,
        help="Clusters list file. Default: <outdir>/clusters.list",
    )
    ap.add_argument(
        "--jac",
        default=None,
        help="Jaccard cutoff label used to infer default pairs filename (e.g. 0.7 -> pairs.jaccard_ge0.7.tsv)",
    )
    ap.add_argument(
        "--map-out",
        default="cluster2super.tsv",
        help="Output mapping filename under outdir",
    )
    ap.add_argument(
        "--summary-out",
        default="super_summary.tsv",
        help="Output summary filename under outdir",
    )
    ap.add_argument(
        "--prefix",
        default="super_",
        help="Supercluster ID prefix (default: super_)",
    )
    ap.add_argument(
        "--pad",
        type=int,
        default=4,
        help="Zero-padding width for supercluster IDs (default: 4 -> super_0001)",
    )
    return ap.parse_args()


args = parse_args()
outdir = Path(args.outdir)

if args.pairs:
    pairs = Path(args.pairs)
else:
    if args.jac is not None:
        pairs = outdir / f"pairs.jaccard_ge{args.jac}.tsv"
    else:
        pairs = outdir / "pairs.jaccard.tsv"

clusters_file = Path(args.clusters) if args.clusters else (outdir / "clusters.list")
clusters = [x.strip() for x in clusters_file.read_text().splitlines() if x.strip()]

parent = {c: c for c in clusters}
rank = {c: 0 for c in clusters}

def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x

def union(a, b):
    ra, rb = find(a), find(b)
    if ra == rb:
        return
    if rank[ra] < rank[rb]:
        parent[ra] = rb
    elif rank[ra] > rank[rb]:
        parent[rb] = ra
    else:
        parent[rb] = ra
        rank[ra] += 1

if pairs.exists():
    for line in pairs.read_text().splitlines():
        if not line.strip():
            continue
        a, b, _ = line.split("\t")
        union(a, b)

groups = {}
for c in clusters:
    r = find(c)
    groups.setdefault(r, []).append(c)

roots = sorted(groups.keys())

map_out = outdir / args.map_out
summary_out = outdir / args.summary_out

with open(map_out, "w") as f:
    for i, r in enumerate(roots, 1):
        sid = f"{args.prefix}{i:0{args.pad}d}"
        for c in sorted(groups[r]):
            f.write(f"{c}\t{sid}\n")

with open(summary_out, "w") as f:
    for i, r in enumerate(roots, 1):
        f.write(
            f"{args.prefix}{i:0{args.pad}d}\t{len(groups[r])}\t{','.join(sorted(groups[r]))}\n"
        )

print(f"OK: wrote {map_out.name} and {summary_out.name}")
