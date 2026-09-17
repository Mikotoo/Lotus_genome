#!/usr/bin/env bash
set -euo pipefail

############################################
# 输入文件（请按需修改）
############################################
COORDS=GifuT2T_MG20T2T.filtered.coords

Gifu_T2T=Lotus_GifuT2T_v1.0.fasta
MG20_T2T=Lotus_MG20T2T_v1.0.fasta

PREFIX=GifuT2T_vs_MG20T2T

############################################
# 0. 准备 genome 文件
############################################
samtools faidx "$Gifu_T2T"
samtools faidx "$MG20_T2T"

cut -f1,2 ${Gifu_T2T}.fai > GifuT2T.genome
cut -f1,2 ${MG20_T2T}.fai > MG20T2T.genome

############################################
# 1. 提取 GifuT2T（query）比对区段
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
> ${PREFIX}.GifuT2T.aligned.bed

############################################
# 2. 提取 MG20T2T（ref）比对区段
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
> ${PREFIX}.MG20T2T.aligned.bed

############################################
# 3. 反选（complement）得到 PAV
############################################
bedtools complement \
  -i ${PREFIX}.GifuT2T.aligned.bed \
  -g GifuT2T.genome \
> ${PREFIX}.GifuT2T.PAV.bed

bedtools complement \
  -i ${PREFIX}.MG20T2T.aligned.bed \
  -g MG20T2T.genome \
> ${PREFIX}.MG20T2T.PAV.bed

############################################
# 4. 统计区间数量 & 总长度
############################################

for f in ${PREFIX}.GifuT2T.PAV.bed ${PREFIX}.MG20T2T.PAV.bed
do

    echo "File: $f" >> ${PREFIX}.PAV.summary

    awk '
    {
        len = $3 - $2 + 1
        if (len > 0) {
            n++
            sum += len
        }
    }
    END {
        printf("Intervals: %d\nTotal length (bp): %d\n", n, sum)
    }' "$f" >> ${PREFIX}.PAV.summary
done
