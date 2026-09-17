#!/usr/bin/env bash
set -euo pipefail

FEATURE_BED="gifuT2T_newseq.bed"
BIN_BED="{PROJ_04_LOTUS_GENOME}/00_assembly_result/GifuT2T/Gifu_100k.bed"
OUT_TSV="gifu_100kb_coverage.tsv"

# 1) 规范化：只取前三列，确保start<end，排序
awk 'BEGIN{OFS="\t"}
     NF>=3 && $2<$3 {print $1,$2,$3}' "$FEATURE_BED" \
  | sort -k1,1 -k2,2n > features.sorted.bed

awk 'BEGIN{OFS="\t"}
     NF>=3 && $2<$3 {print $1,$2,$3}' "$BIN_BED" \
  | sort -k1,1 -k2,2n > bins100k.sorted.bed

# 2) 计算 features 与 bin 的重叠长度（-wo 最后列是 overlap bp）
# 输出列：
# bin_chr bin_start bin_end feat_chr feat_start feat_end overlap
bedtools intersect -a bins100k.sorted.bed -b features.sorted.bed -wo \
  | awk 'BEGIN{OFS="\t"}
         {
           key=$1"\t"$2"\t"$3;
           ov[key]+=$7;
         }
         END{
           for(k in ov){
             print k, ov[k];
           }
         }' \
  | sort -k1,1 -k2,2n > bin_overlap_bp.tsv

# 3) 把没有任何 overlap 的 bin 补 0，并计算比例 = overlap / bin_len
# 输出：chr start end bin_len overlap_bp proportion
awk 'BEGIN{OFS="\t"}
     FNR==NR{
       ov[$1"\t"$2"\t"$3]=$4;
       next
     }
     {
       key=$1"\t"$2"\t"$3;
       len=$3-$2;
       bp=(key in ov)?ov[key]:0;
       prop=(len>0)?bp/len:0;
       print $1,$2,$3,len,bp,prop
     }' bin_overlap_bp.tsv bins100k.sorted.bed > "$OUT_TSV"

echo "Done: $OUT_TSV"
