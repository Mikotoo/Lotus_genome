#!/usr/bin/env python3
"""
Filter predicted proteins and corresponding CDS and GFF3 features:
    - Keep proteins >= min length (default 50 aa), starting with M, and with no internal stop codons.
    - Output filtered proteins and (optionally) matching CDS.
    - Output GFF3 with removed features whose ID or Parent matches filtered-out proteins.
"""

import argparse
from Bio import SeqIO
import re
from collections import defaultdict

def parse_args():
    ap = argparse.ArgumentParser(
        description="Filter proteins by length, N-terminal methionine, and premature stop codon; filter corresponding CDS and GFF3."
    )
    ap.add_argument("-i", "--input", required=True, help="protein fasta")
    ap.add_argument("-o", "--out", required=True, help="output protein fasta")
    ap.add_argument("-l", "--length", type=int, default=20, help="minimum protein length")
    ap.add_argument("-c", "--cds", required=False, help="input CDS fasta (optional)")
    ap.add_argument("-C", "--out_cds", required=False, help="output filtered CDS fasta (optional)")
    ap.add_argument("-g", "--gff", required=False, help="input GFF3 (optional)")
    ap.add_argument("-G", "--out_gff", required=False, help="output filtered GFF3 (optional)")
    return ap.parse_args()

def protein_pass(protein_seq, min_length):
    """Return True if sequence passes all criteria."""
    s = str(protein_seq)
    if len(s) < min_length:
        return False
    if s[0] != "M":
        return False
    # Allow only a single terminal stop, no internal
    if "*" in s[:-1]:
        return False
    return True

def filter_proteins_and_write(input_fasta, min_length, output_fasta):
    """Return set of IDs of proteins that pass filters and write them."""
    passing_ids = set()
    with open(output_fasta, "w") as out:
        for rec in SeqIO.parse(input_fasta, "fasta"):
            if protein_pass(rec.seq, min_length):
                SeqIO.write(rec, out, "fasta")
                passing_ids.add(rec.id.split()[0])
    return passing_ids

def filter_cds_and_write(cds_fasta, passing_ids, output_fasta):
    n = 0
    with open(output_fasta, "w") as out:
        for rec in SeqIO.parse(cds_fasta, "fasta"):
            if rec.id.split()[0] in passing_ids:
                SeqIO.write(rec, out, "fasta")
                n += 1
    return n

def parse_attr(attr_field):
    """Parse GFF3 attribute field into a dict."""
    attrs = {}
    for part in attr_field.strip().split(";"):
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            attrs[k] = v
    return attrs

def collect_transcript_gene_relations(gff_file):
    """Return mappings: transcript → gene and gene → [transcripts]."""
    t2g = {}
    g2t = defaultdict(list)
    with open(gff_file) as fh:
        for ln in fh:
            if ln.startswith("#") or ln.strip() == "":
                continue
            flds = ln.rstrip().split("\t")
            if len(flds) < 9:
                continue
            ftype, attrs = flds[2], flds[8]
            a = parse_attr(attrs)
            if ftype.lower() in {"mrna", "transcript"}:
                tid = a.get("ID")
                gid = a.get("Parent")
                if tid and gid:
                    t2g[tid] = gid
                    g2t[gid].append(tid)
    return t2g, g2t

def decide_gene_removal(g2t, filtered_out_ids):
    """Genes to drop = every transcript filtered out."""
    drop_genes = set()
    for gid, tids in g2t.items():
        if all(tid in filtered_out_ids for tid in tids):
            drop_genes.add(gid)
    return drop_genes

def filter_gff(gff_in, gff_out, passing_ids, drop_genes):
    kept, removed = 0, 0
    with open(gff_in) as inp, open(gff_out, "w") as out:
        for ln in inp:
            if ln.startswith("#") or ln.strip() == "":
                out.write(ln)
                continue
            flds = ln.rstrip().split("\t")
            if len(flds) < 9:
                out.write(ln)
                continue
            ftype, attrs = flds[2], flds[8]
            a = parse_attr(attrs)
            gid = a.get("ID") if ftype.lower() == "gene" else None
            tid = None
            if ftype.lower() in {"mrna", "transcript"}:
                tid = a.get("ID")
            elif "Parent" in a:
                tid = a["Parent"].split(",")[0]  # first parent

            # drop rules
            # Remove gene if all transcripts filtered out
            if gid and gid in drop_genes:
                removed += 1
                continue
            # Remove transcript and all features whose Parent matches
            if tid and tid not in passing_ids:
                removed += 1
                continue

            out.write(ln)
            kept += 1
    return kept, removed

def main():
    args = parse_args()
    passing_ids = filter_proteins_and_write(args.input, args.length, args.out)
    print(f"{len(passing_ids)} proteins passed all filters and written to {args.out}")
    if args.cds and args.out_cds:
        n = filter_cds_and_write(args.cds, passing_ids, args.out_cds)
        print(f"{n} CDS records matched filtered proteins and written to {args.out_cds}")
    if args.gff and args.out_gff:
        t2g, g2t = collect_transcript_gene_relations(args.gff)
        filtered_out_ids = set(t2g.keys()) - passing_ids
        drop_genes = decide_gene_removal(g2t, filtered_out_ids)
        kept, removed = filter_gff(args.gff, args.out_gff, passing_ids, drop_genes)
        print(f"GFF filtered: kept {kept} lines; removed {removed} lines → {args.out_gff}")

if __name__ == "__main__":
    main()
