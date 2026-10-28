#!/usr/bin/env bash
set -euo pipefail

############################################
# 输入文件（请按需修改）
############################################
COORDS=GifuT2T_GifuOld.filtered.coords

FA_T2T=Lotus_GifuT2T_v1.0.fasta
FA_OLD=Lotus_Gifu_Old1.fa

PREFIX=GifuT2T_vs_GifuOld

############################################
# 0. 准备 genome 文件
############################################
samtools faidx "$FA_T2T"
samtools faidx "$FA_OLD"

cut -f1,2 ${FA_T2T}.fai > GifuT2T.genome
cut -f1,2 ${FA_OLD}.fai > GifuOld.genome

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
# 2. 提取 GifuOld（ref）比对区段
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
> ${PREFIX}.GifuOld.aligned.bed

############################################
# 3. 反选（complement）得到 PAV
############################################
bedtools complement \
  -i ${PREFIX}.GifuT2T.aligned.bed \
  -g GifuT2T.genome \
> ${PREFIX}.GifuT2T.PAV.bed

bedtools complement \
  -i ${PREFIX}.GifuOld.aligned.bed \
  -g GifuOld.genome \
> ${PREFIX}.GifuOld.PAV.bed

############################################
# 4. 统计区间数量 & 总长度
############################################

for f in ${PREFIX}.GifuT2T.PAV.bed ${PREFIX}.GifuOld.PAV.bed
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
