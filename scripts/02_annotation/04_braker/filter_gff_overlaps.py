#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Filter overlapping/duplicate genes from a GFF3, while preserving *all* child features
(mRNA/exon/CDS/UTR/...) of the kept genes and keeping output format/order close to input.

Usage:
  python filter_gff_overlaps_full.py in.gff3 out.gff3 \
    --strand-aware \
    --reciprocal-threshold 0.0 \
    --keep-criterion longest

Options:
  --strand-aware            Only treat same-strand genes as overlapping (default off).
  --reciprocal-threshold    Reciprocal overlap threshold in [0,1]; 0 = any overlap (default 0).
  --keep-criterion          Which gene to keep within a cluster: longest | id-min (default longest).
"""

import sys
import argparse
from collections import defaultdict, OrderedDict, deque

# -------------------- GFF utils --------------------
def parse_attrs(attr_field: str) -> OrderedDict:
    attrs = OrderedDict()
    for entry in attr_field.split(";"):
        entry = entry.strip()
        if not entry:
            continue
        if "=" in entry:
            k, v = entry.split("=", 1)
            attrs[k] = v
        else:
            attrs[entry] = ""
    return attrs

def attrs_to_str(attrs: OrderedDict) -> str:
    return ";".join(f"{k}={v}" if v != "" else k for k, v in attrs.items())

def length1(s: int, e: int) -> int:
    return e - s + 1  # GFF is 1-based inclusive

def overlap_len(a, b) -> int:
    # a,b are dicts with start,end
    return max(0, min(a["end"], b["end"]) - max(a["start"], b["start"]) + 1)

def reciprocal_overlap(a, b) -> float:
    inter = overlap_len(a, b)
    if inter == 0:
        return 0.0
    la = length1(a["start"], a["end"])
    lb = length1(b["start"], b["end"])
    return min(inter / la, inter / lb)

# -------------------- Parsing whole GFF --------------------
def read_gff(path):
    """
    Return:
      headers: list[str]                 # lines starting with '#' or malformed kept as-is
      records: dict[id]->rec             # every feature with ID (gene/mRNA/exon/CDS/... if they have ID)
      anon_lines: list[(lineno,line)]    # features without ID (会作为“匿名记录”保留，但仅在其 Parent 最终被保留时输出)
      parent_of: dict[id]->list[parent_ids]
      children_of: dict[id]->list[child_ids]
      genes: dict[gid]->rec              # only feature type == 'gene'
      order_index: dict[id]->lineno      # first appearance order for features with ID
    """
    headers = []
    records = {}          # id -> rec
    anon_lines = []       # lines without ID but may have Parent (exon 等常有 ID，也可能没有)
    parent_of = defaultdict(list)    # child_id -> [parent_id,...]
    children_of = defaultdict(list)  # parent_id -> [child_id,...]
    genes = {}
    order_index = {}

    with open(path, "r") as fh:
        lineno = 0
        for raw in fh:
            lineno += 1
            line = raw.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                headers.append(line)
                continue

            parts = line.split("\t")
            if len(parts) < 9:
                # 非标准行，当作 header 原样保留
                headers.append(line)
                continue

            seqid, source, ftype, start, end, score, strand, phase, attr = parts
            try:
                start_i = int(start); end_i = int(end)
            except Exception:
                headers.append(line)
                continue

            attrs = parse_attrs(attr)
            fid = attrs.get("ID", None)
            parents = []
            if "Parent" in attrs:
                parents = [p for p in attrs["Parent"].split(",") if p]

            rec = {
                "chr": seqid,
                "source": source,
                "type": ftype,
                "start": start_i,
                "end": end_i,
                "score": score,
                "strand": strand if strand in "+-" else ".",
                "phase": phase,
                "attrs": attrs,
                "line": line,
                "lineno": lineno,
            }

            if fid:
                records[fid] = rec
                if fid not in order_index:
                    order_index[fid] = lineno
                for p in parents:
                    parent_of[fid].append(p)
                    children_of[p].append(fid)
                if ftype.lower() == "gene":
                    genes[fid] = rec
            else:
                # 没有 ID 的行：先收集，等到输出阶段根据其 Parent 是否被保留决定是否写出
                anon_lines.append((lineno, rec, parents))

    return headers, records, anon_lines, parent_of, children_of, genes, order_index

# -------------------- Dedup & clustering on genes --------------------
def dedup_exact(genes_items):
    """同 chr,start,end,strand 完全一致视为重复，只保留一个（按最长、再按 ID 排序）。"""
    buckets = defaultdict(list)
    for gid, g in genes_items:
        key = (g["chr"], g["start"], g["end"], g["strand"])
        buckets[key].append((gid, g))

    kept, dropped = [], set()
    for key, lst in buckets.items():
        if len(lst) == 1:
            kept.append(lst[0])
        else:
            lst.sort(key=lambda x: (-(x[1]["end"] - x[1]["start"]), x[0]))
            kept.append(lst[0])
            for gid, _ in lst[1:]:
                dropped.add(gid)
    return kept, dropped

def build_gene_clusters(genes_items, strand_aware=False, reciprocal_threshold=0.0):
    """把基因连通为簇。"""
    groups = defaultdict(list)  # per (chr[,strand])
    for gid, g in genes_items:
        key = (g["chr"], g["strand"]) if strand_aware else (g["chr"], None)
        groups[key].append((gid, g))

    clusters = []
    for key, arr in groups.items():
        arr.sort(key=lambda x: x[1]["start"])
        n = len(arr)
        visited = [False]*n
        for i in range(n):
            if visited[i]:
                continue
            comp = []
            dq = deque([i]); visited[i] = True
            while dq:
                u = dq.popleft()
                comp.append(u)
                urec = arr[u][1]
                # 只需检查 start <= u.end 的右侧邻域
                v = u + 1
                while v < n and arr[v][1]["start"] <= urec["end"]:
                    if not visited[v]:
                        vrec = arr[v][1]
                        ok = False
                        if reciprocal_threshold <= 0:
                            ok = not (vrec["start"] > urec["end"] or vrec["end"] < urec["start"])
                        else:
                            ok = (reciprocal_overlap(urec, vrec) >= reciprocal_threshold)
                        if ok:
                            visited[v] = True
                            dq.append(v)
                    v += 1
            clusters.append([arr[idx] for idx in comp])
    return clusters

def choose_rep(cluster, keep_criterion="longest"):
    if keep_criterion == "id-min":
        cluster.sort(key=lambda x: x[0])
        return cluster[0]
    # longest
    cluster.sort(key=lambda x: (-(x[1]["end"] - x[1]["start"] + 1), x[0]))
    return cluster[0]

# -------------------- Collect descendants & rewrite --------------------
def collect_descendants(root_gene_id, children_of, records):
    """返回以 gene 为根的所有后代（包含 gene 自身）的 ID 集合。"""
    keep_ids = set()
    stack = [root_gene_id]
    while stack:
        pid = stack.pop()
        if pid in keep_ids:
            continue
        keep_ids.add(pid)
        for cid in children_of.get(pid, []):
            if cid in records:       # 只收有 ID 的
                stack.append(cid)
    return keep_ids

def prune_parent_attr(rec, kept_ids_set):
    """修剪 Parent= 使之只保留仍存在的父 ID；若清空则返回 False 表示应丢弃该行。"""
    attrs = rec["attrs"]
    if "Parent" not in attrs:
        return True
    parents = [p for p in attrs["Parent"].split(",") if p]
    new_parents = [p for p in parents if p in kept_ids_set]
    if not new_parents:
        return False
    if ",".join(new_parents) != attrs["Parent"]:
        attrs = attrs.copy()
        attrs["Parent"] = ",".join(new_parents)
        rec["attrs"] = attrs
    return True

def rec_to_line(rec):
    # 用当前 attrs 重建第 9 列，并保持其他列不变
    parts = [
        rec["chr"],
        rec["source"],
        rec["type"],
        str(rec["start"]),
        str(rec["end"]),
        rec["score"],
        rec["strand"],
        rec["phase"],
        attrs_to_str(rec["attrs"])
    ]
    return "\t".join(parts)

# -------------------- Main --------------------
def main():
    ap = argparse.ArgumentParser(description="Filter overlapping genes from GFF3 and keep full structures.")
    ap.add_argument("input")
    ap.add_argument("output")
    ap.add_argument("--strand-aware", action="store_true", help="Overlap only within same strand.")
    ap.add_argument("--reciprocal-threshold", type=float, default=0.0,
                    help="Reciprocal overlap threshold in [0,1]. 0: any overlap.")
    ap.add_argument("--keep-criterion", choices=["longest","id-min"], default="longest",
                    help="Representative in a cluster.")
    args = ap.parse_args()

    headers, records, anon_lines, parent_of, children_of, genes, order_index = read_gff(args.input)

    # 1) 去完全重复
    gene_items = [(gid, genes[gid]) for gid in genes]
    dedup_kept, dropped_exact = dedup_exact(gene_items)

    # 2) 基于重叠建立簇并选择代表
    clusters = build_gene_clusters(dedup_kept,
                                   strand_aware=args.strand_aware,
                                   reciprocal_threshold=args.reciprocal_threshold)
    kept_gene_ids = set()
    for cl in clusters:
        gid, _ = choose_rep(cl, keep_criterion=args.keep_criterion)
        kept_gene_ids.add(gid)

    # 3) 收集每个保留基因的所有后代（递归）
    kept_ids = set()
    for gid in kept_gene_ids:
        kept_ids |= collect_descendants(gid, children_of, records)
    # 顺带把基因自身也包含在 kept_ids 里（collect 已含）

    # 4) 对匿名行（无 ID，但有 Parent 的），仅当 Parent 之一被保留时写出，并修剪 Parent
    #    对有 ID 的行：如果在 kept_ids，则写出；否则丢弃
    output_items = []  # (lineno, line_str)
    # 4.1 有 ID 的
    for fid in kept_ids:
        if fid not in records:
            continue
        rec = records[fid]
        # 修剪 Parent（例如 mRNA 的 Parent=gene；若 gene 被保留则 OK；若还有其他被删的 Parent 则去掉）
        if not prune_parent_attr(rec, kept_ids):
            continue
        # 只对 Parent 变动的行重构；否则尽量用原行，减少改动
        line_out = rec_to_line(rec) if "Parent" in rec["attrs"] and line_parent_changed(rec["line"], rec["attrs"]) else rec["line"]
        output_items.append((rec["lineno"], line_out))

    # 4.2 无 ID 的（多数情况下也会有 Parent）
    for lineno, rec, parents in anon_lines:
        if not parents:
            # 无 Parent 的匿名行：通常与结构无关，保留原样
            output_items.append((lineno, rec["line"]))
        else:
            keep = [p for p in parents if p in kept_ids]
            if keep:
                # 修剪 Parent
                if "attrs" in rec and "Parent" in rec["attrs"]:
                    attrs = rec["attrs"].copy()
                    attrs["Parent"] = ",".join(keep)
                    rec2 = rec.copy(); rec2["attrs"] = attrs
                    output_items.append((lineno, rec_to_line(rec2)))
                else:
                    output_items.append((lineno, rec["line"]))

    # 5) 拼装输出（头注释 + 按原行号排序后的主体；若存在 ##FASTA 段在 headers 里也会被保留）
    output_items.sort(key=lambda x: x[0])
    with open(args.output, "w") as out:
        for h in headers:
            out.write(h + "\n")
        for _, ln in output_items:
            out.write(ln + "\n")

    # 统计
    sys.stderr.write(f"[INFO] Input genes: {len(genes)}\n")
    sys.stderr.write(f"[INFO] Removed exact duplicates: {len(dropped_exact)}\n")
    sys.stderr.write(f"[INFO] Kept genes (after overlap filtering): {len(kept_gene_ids)}\n")
    sys.stderr.write(f"[INFO] Total features written: {len(output_items)}\n")
    sys.stderr.write(f"[INFO] Output: {args.output}\n")

def line_parent_changed(orig_line: str, new_attrs: OrderedDict) -> bool:
    """轻量判断 Parent 是否被修改，决定是否重建第9列（尽量减少改动）。"""
    try:
        old_attr = orig_line.split("\t", 8)[8]
    except Exception:
        return True
    old = parse_attrs(old_attr)
    return old.get("Parent","") != new_attrs.get("Parent","")

if __name__ == "__main__":
    main()

