#!/bin/bash

# 设置文件路径
new_seq_bed="gifuT2T_newseq.bed"
cent_bed="cent.bed"
rDNA_bed="rDNA.bed"
tel_bed="gifu_tel.bed"
TE_bed="repeat2.bed"

# 计算新序列与每个区域的重叠部分（去重后的TE区域）
cent_overlap=$(bedtools intersect -a $new_seq_bed -b $cent_bed -wa | awk '{sum+=$3-$2} END {print sum}')
rDNA_overlap=$(bedtools intersect -a $new_seq_bed -b $rDNA_bed -wa | awk '{sum+=$3-$2} END {print sum}')
tel_overlap=$(bedtools intersect -a $new_seq_bed -b $tel_bed -wa | awk '{sum+=$3-$2} END {print sum}')

# 计算新序列与TE的重叠区域并去重
TE_overlap_raw=$(bedtools intersect -a $new_seq_bed -b $TE_bed -wa)
TE_overlap_unique=$(echo "$TE_overlap_raw" | bedtools sort | bedtools merge -i - | awk '{sum+=$3-$2} END {print sum}')

# 计算新序列的总长度
new_seq_length=$(awk '{sum+=$3-$2} END {print sum}' $new_seq_bed)

# 计算其他区域的长度（不在cent, rDNA, tel, TE内的区域）
other_overlap=$(($new_seq_length - $cent_overlap - $rDNA_overlap - $tel_overlap - $TE_overlap_unique))

# 输出占比
echo "Centromere overlap: $cent_overlap bp"
echo "rDNA overlap: $rDNA_overlap bp"
echo "Telomere overlap: $tel_overlap bp"
echo "TE overlap (unique): $TE_overlap_unique bp"
echo "Other regions overlap: $other_overlap bp"

# 使用Python绘制饼图并保存为PDF文件
python3 - <<EOF
import matplotlib.pyplot as plt

# 占比数据
labels = ['Centromere', 'rDNA', 'Telomere', 'Transposable Elements', 'Other Regions']
sizes = [$cent_overlap, $rDNA_overlap, $tel_overlap, $TE_overlap_unique, $other_overlap]
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
