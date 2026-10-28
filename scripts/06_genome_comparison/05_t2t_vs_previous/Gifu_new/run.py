import matplotlib.pyplot as plt

def read_bed(file_path):
    """读取BED文件，返回列表[(chrom, start, end)]"""
    regions = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('#') or line.strip() == "":
                continue
            chrom, start, end = line.strip().split()[:3]
            regions.append((chrom, int(start), int(end)))
    return regions

def calculate_overlap(region1, region2):
    """计算两个区域的重叠部分长度"""
    chrom1, start1, end1 = region1
    chrom2, start2, end2 = region2
    # 确保是同一个染色体
    if chrom1 != chrom2:
        return 0
    # 计算重叠区间
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    if overlap_start < overlap_end:
        return overlap_end - overlap_start
    return 0

def subtract_region(region1, region2):
    """从region1中去除region2的部分"""
    chrom1, start1, end1 = region1
    chrom2, start2, end2 = region2
    if chrom1 != chrom2:
        return [region1]
    result = []
    # 如果region1完全在region2之外，则直接返回
    if end1 <= start2 or start1 >= end2:
        return [region1]
    # 计算去除后的区域
    if start1 < start2:
        result.append((chrom1, start1, start2))
    if end1 > end2:
        result.append((chrom1, end2, end1))
    return result

def calculate_total_length(regions):
    """计算区域总长度"""
    return sum([end - start for _, start, end in regions])

# 读取文件
new_seq_bed = "gifuT2T_newseq.bed"
cent_bed = "cent.bed"
rDNA_bed = "rDNA.bed"
tel_bed = "gifu_tel.bed"
TE_bed = "repeat2.bed"

new_seq_regions = read_bed(new_seq_bed)
cent_regions = read_bed(cent_bed)
rDNA_regions = read_bed(rDNA_bed)
tel_regions = read_bed(tel_bed)
TE_regions = read_bed(TE_bed)

# 计算新序列与各个区域的重叠部分的长度
cent_overlap = sum([calculate_overlap(seq, cent) for seq in new_seq_regions for cent in cent_regions])
rDNA_overlap = sum([calculate_overlap(seq, rDNA) for seq in new_seq_regions for rDNA in rDNA_regions])
tel_overlap = sum([calculate_overlap(seq, tel) for seq in new_seq_regions for tel in tel_regions])
TE_overlap = sum([calculate_overlap(seq, TE) for seq in new_seq_regions for TE in TE_regions])

# 计算新序列的总长度
new_seq_length = calculate_total_length(new_seq_regions)

# 逐步去除每个区域的交集部分，并计算剩余的部分
remaining_after_cent = new_seq_regions
for cent in cent_regions:
    remaining_after_cent = [seq for seq in remaining_after_cent if calculate_overlap(seq, cent) == 0]

remaining_after_rDNA = remaining_after_cent
for rDNA in rDNA_regions:
    remaining_after_rDNA = [seq for seq in remaining_after_rDNA if calculate_overlap(seq, rDNA) == 0]

remaining_after_tel = remaining_after_rDNA
for tel in tel_regions:
    remaining_after_tel = [seq for seq in remaining_after_tel if calculate_overlap(seq, tel) == 0]

remaining_after_TE = remaining_after_tel
for TE in TE_regions:
    remaining_after_TE = [seq for seq in remaining_after_TE if calculate_overlap(seq, TE) == 0]

# 计算其他区域的长度
other_overlap = calculate_total_length(remaining_after_TE)

# 计算占比
cent_ratio = cent_overlap / new_seq_length
rDNA_ratio = rDNA_overlap / new_seq_length
tel_ratio = tel_overlap / new_seq_length
TE_ratio = TE_overlap / new_seq_length
other_ratio = other_overlap / new_seq_length

# 输出各个区域的交集和占比
print(f"Centromere overlap: {cent_overlap} bp, Centromere ratio: {cent_ratio:.4f}")
print(f"rDNA overlap: {rDNA_overlap} bp, rDNA ratio: {rDNA_ratio:.4f}")
print(f"Telomere overlap: {tel_overlap} bp, Telomere ratio: {tel_ratio:.4f}")
print(f"TE overlap (unique): {TE_overlap} bp, TE ratio: {TE_ratio:.4f}")
print(f"Other regions overlap: {other_overlap} bp, Other regions ratio: {other_ratio:.4f}")

# 使用Python绘制饼图并保存为PDF文件
labels = ['Centromere', 'rDNA', 'Telomere', 'Transposable Elements', 'Other Regions']
sizes = [cent_overlap, rDNA_overlap, tel_overlap, TE_overlap, other_overlap]
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
