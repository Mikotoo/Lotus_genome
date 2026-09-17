#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
from xml.etree import ElementTree as ET
import re

def parse_args():
    parser = argparse.ArgumentParser(description="从 BLAST XML 提取物种名称并统计比例（保留前十个命中）")
    parser.add_argument("-i", "--input", required=True, help="输入的 BLAST XML 文件")
    parser.add_argument("-a", "--accession2taxid", required=True, help="accession2taxid 文件")
    parser.add_argument("-n", "--names", required=True, help="names.dmp 文件")
    parser.add_argument("-o", "--output_prefix", required=True, help="输出前缀")
    return parser.parse_args()

def extract_top_hits(xml_file, output_tiqugi_file):
    tree = ET.parse(xml_file)
    root = tree.getroot()
    hit_list = []
    
    with open(output_tiqugi_file, "w") as out:
        for iteration in root.findall(".//Iteration"):
            query_id = iteration.findtext("Iteration_query-def")
            hits = iteration.findall("Iteration_hits/Hit")
            for hit in hits[:1]:  # 限制最多取前十个
                hit_id = hit.findtext("Hit_id") or "-"
                hit_def = hit.findtext("Hit_def") or "-"
                hit_list.append((query_id, hit_id, hit_def))
                out.write(f"{query_id}\t{hit_id}\t{hit_def}\n")
    
    return hit_list

def load_accession2taxid(file):
    acc2tax = {}
    with open(file) as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split("\t")  # 确保按制表符分割字段
            if len(parts) >= 2:
                taxid, accession = parts[0], parts[1]
                # 去除 accession 的前缀（例如 ref|, gb| 等）
                accession = re.sub(r"^(ref|gb)\|", "", accession)
                acc2tax[accession] = taxid
    return acc2tax

def load_taxid2name(file):
    tax2name = {}
    with open(file) as f:
        for line in f:
            parts = [p.strip() for p in line.strip().split("|")]
            if len(parts) >= 4 and parts[3] == "scientific name":
                taxid = parts[0]
                name = parts[1]
                tax2name[taxid] = name
    return tax2name

def annotate_and_write(hit_list, acc2tax, tax2name, output_annotated_file):
    total = 0
    matched = 0
    species_counter = {}
    
    with open(output_annotated_file, "w") as out:
        for query_id, hit_id, hit_def in hit_list:
            parts = hit_id.split("|")
            accession = parts[1] if len(parts) > 1 else None
            if accession and accession in acc2tax:
                taxid = acc2tax[accession]
                species = tax2name.get(taxid)
                if species:
                    matched += 1
                    species_counter[species] = species_counter.get(species, 0) + 1
                    out.write(f"{query_id}\t{hit_id}\t{hit_def}\t{species}\n")
            total += 1
    return species_counter, total

def write_summary(species_counter, total, output_summary_file):
    with open(output_summary_file, "w") as out:
        out.write("Name\tHit_reads\tpercentage1\tpercentage2\n")
        for species, count in species_counter.items():
            perc1 = f"{(count / total) * 100:.2f}%"
            perc2 = f"{count / total:.4f}"
            out.write(f"{species}\t{count}\t{perc1}\t{perc2}\n")

def main():
    args = parse_args()

    print("📥 提取 BLAST XML 中的命中信息 ...")
    tiqugi_file = args.output_prefix + "_best_tiqu_gi.txt"
    hits = extract_top_hits(args.input, tiqugi_file)

    print("📄 加载 accession → taxid 映射 ...")
    acc2tax = load_accession2taxid(args.accession2taxid)

    print("📄 加载 taxid → scientific name 映射 ...")
    tax2name = load_taxid2name(args.names)

    print("🧬 注释物种名并写出 ...")
    annotated_file = args.output_prefix + "_scientific_name.txt"
    species_counter, total = annotate_and_write(hits, acc2tax, tax2name, annotated_file)

    print("📊 统计物种比例并写出 ...")
    summary_file = args.output_prefix + "_summary.txt"
    write_summary(species_counter, total, summary_file)

    print(f"✅ 完成！输出文件：\n{tiqugi_file}\n{annotated_file}\n{summary_file}")

if __name__ == "__main__":
    main()
