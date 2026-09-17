#!/usr/bin/env bash
set -euo pipefail

############################################
# 输入文件（请按需修改）
############################################
COORDS=MG20T2T_MG20Old.filtered.coords

FA_T2T=Lotus_MG20T2T_v1.0.fasta
FA_OLD=Lotus_MG20_Old1.fa

PREFIX=MG20T2T_vs_MG20Old

############################################
# 0. 准备 genome 文件
############################################
samtools faidx "$FA_T2T"
samtools faidx "$FA_OLD"

cut -f1,2 ${FA_T2T}.fai > MG20T2T.genome
cut -f1,2 ${FA_OLD}.fai > MG20Old.genome

############################################
# 1. 提取 MG20T2T（query）比对区段
############################################
awk -F'\t' '
{
    chr = $11
    s = $3
    e = $4
    if (s > e) {tmp=s; s=e; e=tmp}
    if (e > s)
        print chr "\t" s-1 "\t" e
}
' "$COORDS" \
| sort -k1,1 -k2,2n \
| bedtools merge -d 50 \
> ${PREFIX}.MG20T2T.aligned.bed

############################################
# 2. 提取 MG20Old（ref）比对区段
############################################
awk -F'\t' '
{
    chr = $10
    s = $1
    e = $2
    if (s > e) {tmp=s; s=e; e=tmp}
    if (e > s)
        print chr "\t" s-1 "\t" e
}
' "$COORDS" \
| sort -k1,1 -k2,2n \
| bedtools merge -d 50 \
> ${PREFIX}.MG20Old.aligned.bed

############################################
# 3. 反选（complement）得到 PAV
############################################
bedtools complement \
  -i ${PREFIX}.MG20T2T.aligned.bed \
  -g MG20T2T.genome \
> ${PREFIX}.MG20T2T.PAV.bed

bedtools complement \
  -i ${PREFIX}.MG20Old.aligned.bed \
  -g MG20Old.genome \
> ${PREFIX}.MG20Old.PAV.bed

############################################
# 4. 统计区间数量 & 总长度
############################################

for f in ${PREFIX}.MG20T2T.PAV.bed ${PREFIX}.MG20Old.PAV.bed
do

    echo "File: $f" >> ${PREFIX}.PAV.summary

    awk '
    {
        len = $3 - $2
        if (len > 0) {
            n++
            sum += len
        }
    }
    END {
        printf("Intervals: %d\nTotal length (bp): %d\n", n, sum)
    }' "$f" >> ${PREFIX}.PAV.summary
done
