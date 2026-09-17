#!/bin/bash
#CSUB -J variation
#CSUB -q c01
#CSUB -o variation.out
#CSUB -e variation.error
#CSUB -n 64
#CSUB -R span[hosts=1]

# 0. 去除原本基因组中的多余的染色体
seqkit grep -v -r -p "LjG1.1_chloro" -p "LjG1.1_chr0" -p "LjG1.1_mito" \
 Lotus_Gifu_Old.fa > Lotus_Gifu_Old1.fa

# 1. 全基因组比对（根据需要调整参数）
# 使用nucmer进行比对，参数可根据前述讨论选择 --mum 或 --maxmatch
# 这里示例采用更全面的--maxmatch并设置较高阈值避免噪音
nucmer -c 1000 --mum --maxgap=1000 --prefix=GifuT2T_GifuOld Lotus_Gifu_Old1.fa Lotus_GifuT2T_v1.0.fasta

# 2. 过滤比对结果，提高比对质量
delta-filter -i 90 -l 1000 -q GifuT2T_GifuOld.delta > GifuT2T_GifuOld.filtered.delta

# 可选：输出比对统计和绘图（非必要步骤，仅检查比对情况）
dnadiff -d GifuT2T_GifuOld.filtered.delta -p GifuT2T_GifuOld.filtered
# show-coords 可用于生成比对坐标表格
show-coords -THrd GifuT2T_GifuOld.filtered.delta > GifuT2T_GifuOld.filtered.coords

# 3. 使用SyRI识别变异
syri -c GifuT2T_GifuOld.filtered.coords -d GifuT2T_GifuOld.filtered.delta \
    -r Lotus_Gifu_Old1.fa -q Lotus_GifuT2T_v1.0.fasta

# 4. 按类别提取变异坐标并保存为BED文件
# 小变异类别：SNP、插入、缺失
awk '$11=="SNP" {print $1"\t"$2"\t"$3 >> "Old_SNP.bed"; print $6"\t"$7"\t"$8 >> "T2T_SNP.bed"}' syri.out
awk '$11=="INS" {print $1"\t"$2"\t"$3 >> "Old_INS.bed"; print $6"\t"$7"\t"$8 >> "T2T_INS.bed"}' syri.out
awk '$11=="DEL" {print $1"\t"$2"\t"$3 >> "Old_DEL.bed"; print $6"\t"$7"\t"$8 >> "T2T_DEL.bed"}' syri.out

# 结构变异类别：倒位、易位、倒置易位、重复、倒置重复
awk '$11=="INV" {print $1"\t"$2"\t"$3 >> "Old_INV.bed"; print $6"\t"$7"\t"$8 >> "T2T_INV.bed"}' syri.out
awk '$11=="TRANS" {print $1"\t"$2"\t"$3 >> "Old_TRANS.bed"; print $6"\t"$7"\t"$8 >> "T2T_TRANS.bed"}' syri.out
awk '$11=="INVTR" {print $1"\t"$2"\t"$3 >> "Old_INVTR.bed"; print $6"\t"$7"\t"$8 >> "T2T_INVTR.bed"}' syri.out
awk '$11=="DUP" {print $1"\t"$2"\t"$3 >> "Old_DUP.bed"; print $6"\t"$7"\t"$8 >> "T2T_DUP.bed"}' syri.out
awk '$11=="INVDP" {print $1"\t"$2"\t"$3 >> "Old_INVDP.bed"; print $6"\t"$7"\t"$8 >> "T2T_INVDP.bed"}' syri.out

# 存在/缺失变异（PAV）：参考或查询中特有的大片段（INS/DEL/NOTAL >50bp）
awk '($11=="INS" || $11=="DEL" || $11=="NOTAL") && (($3-$2+1>50) || ($8-$7+1>50)) \
    {print $1"\t"$2"\t"$3 >> "Old_PAV.bed"; print $6"\t"$7"\t"$8 >> "T2T_PAV.bed"}' syri.out

# 5. 统计各类别变异事件数量和总长度
echo -e "Category\tCount\tTotal_length_in_Old\tTotal_length_in_T2T"
for type in SNP INS DEL INV TRANS INVTR DUP INVDP PAV; do
    # 统计事件数（行数）
    count=$(wc -l < "Old_${type}.bed")
    # 计算参考和查询总长度（bp）
    total_ref=$(awk '{sum += $3-$2+1} END{print sum}' "Old_${type}.bed")
    total_qry=$(awk '{sum += $3-$2+1} END{print sum}' "T2T_${type}.bed")
    echo -e "${type}\t${count}\t${total_ref}\t${total_qry}"
done

#基于syri的结构变异结果
#plotsr --sr syri.out --genomes genomes.txt -o GifuT2T_GifuOld_findSV.pdf

#-s 250000 过滤掉 <250 kb 的小型 rearrangements
#-R 会把没有被结构重排打断的相邻 syntenic blocks 合并起来
plotsr \
    --sr syri.out \
    --genomes genomes.txt \
    -s 250000 \
    -R \
    -W 5.3 \
    -H 5 \
    -f 8 \
    -o GifuT2T_GifuOld_findSV.pdf