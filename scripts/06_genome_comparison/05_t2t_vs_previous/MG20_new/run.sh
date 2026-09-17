#!/bin/bash

# 设置文件路径
new_seq_bed="mg20T2T_newseq.bed"
cent_bed="cent.bed"
rDNA_bed="rDNA.bed"
tel_bed="tel.bed"
TE_bed="repeat2.bed"

# 检查文件路径是否正确
if [[ ! -f "$new_seq_bed" || ! -f "$cent_bed" || ! -f "$rDNA_bed" || ! -f "$tel_bed" || ! -f "$TE_bed" ]]; then
  echo "Error: One or more input files do not exist."
  exit 1
fi

# 计算新序列与cent的交集的真实长度
cent_overlap=$(bedtools intersect -a $new_seq_bed -b $cent_bed -wb | awk '{sum+=$3-$2} END {print sum}')
new_seq_length=$(awk '{sum+=$3-$2} END {print sum}' $new_seq_bed)
cent_ratio=$(echo "scale=4; $cent_overlap / $new_seq_length" | bc)

# 去除cent部分后，计算剩余新序列与rDNA的交集的真实长度
remaining_after_cent=$(bedtools subtract -a $new_seq_bed -b $cent_bed)
rDNA_overlap=$(bedtools intersect -a $remaining_after_cent -b $rDNA_bed -wb | awk '{sum+=$3-$2} END {print sum}')
# 如果rDNA没有交集，则设置为0
if [ -z "$rDNA_overlap" ]; then
    rDNA_overlap=0
fi
rDNA_ratio=$(echo "scale=4; $rDNA_overlap / $new_seq_length" | bc)

# 去除rDNA部分后，计算剩余新序列与telomere的交集的真实长度
remaining_after_rDNA=$(bedtools subtract -a $remaining_after_cent -b $rDNA_bed)
tel_overlap=$(bedtools intersect -a $remaining_after_rDNA -b $tel_bed -wb | awk '{sum+=$3-$2} END {print sum}')
# 如果tel没有交集，则设置为0
if [ -z "$tel_overlap" ]; then
    tel_overlap=0
fi
tel_ratio=$(echo "scale=4; $tel_overlap / $new_seq_length" | bc)

# 去除tel部分后，计算剩余新序列与TE的交集的真实长度
remaining_after_tel=$(bedtools subtract -a $remaining_after_rDNA -b $tel_bed)
TE_overlap_raw=$(bedtools intersect -a $remaining_after_tel -b $TE_bed -wb)
TE_overlap_unique=$(echo "$TE_overlap_raw" | bedtools sort | bedtools merge -i - | awk '{sum+=$3-$2} END {print sum}')
# 如果TE没有交集，则设置为0
if [ -z "$TE_overlap_unique" ]; then
    TE_overlap_unique=0
fi
TE_ratio=$(echo "scale=4; $TE_overlap_unique / $new_seq_length" | bc)

# 最终去除TE后的区域即为Other regions
remaining_after_TE=$(bedtools subtract -a $remaining_after_tel -b $TE_bed)
other_overlap=$(echo "$remaining_after_TE" | awk '{sum+=$3-$2} END {print sum}')
# 如果Other regions没有交集，则设置为0
if [ -z "$other_overlap" ]; then
    other_overlap=0
fi

# 输出各个区域的交集和占比
echo "Centromere overlap: $cent_overlap bp, Centromere ratio: $cent_ratio"
echo "rDNA overlap: $rDNA_overlap bp, rDNA ratio: $rDNA_ratio"
echo "Telomere overlap: $tel_overlap bp, Telomere ratio: $tel_ratio"
echo "TE overlap (unique): $TE_overlap_unique bp, TE ratio: $TE_ratio"
echo "Other regions overlap: $other_overlap bp"

# 使用Python绘制饼图并保存为PDF文件
python3 - <<EOF
import matplotlib.pyplot as plt

# 占比数据
# 如果有空值，将其设置为0
sizes = [$cent_overlap, $rDNA_overlap, $tel_overlap, $TE_overlap_unique, $other_overlap]
labels = ['Centromere', 'rDNA', 'Telomere', 'Transposable Elements', 'Other Regions']
colors = ['#ff9999','#66b3ff','#99ff99','#ffcc99','#c2c2f0']

# 绘制饼图
plt.figure(figsize=(7, 7))
plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=140, colors=colors)
plt.axis('equal')  # 保证饼图是圆形的
plt.title('Proportion of New Sequences in Different Regions')

# 保存为PDF文件
plt.savefig('new_sequences_region_proportions.pdf', format='pdf')

# 显示图形
plt.show()
EOF
