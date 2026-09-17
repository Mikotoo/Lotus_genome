#!/bin/bash
set -euo pipefail

#=================== 参数设置 ===================
cul=$(pwd)
THREADS=88

#=================== 输入文件 ===================
pep1=$1
genome1=$2
gff1=$3
pep2=$4
genome2=$5
gff2=$6
gap=$7   # gap 可为空，但必须占位

#=================== 获取蛋白质序列一致性或覆盖率较低的基因列表 ===================
mkdir -p blastp

# 构建 BLASTP 数据库
makeblastdb -in "$pep1" -dbtype prot -out "blastp/db"

# 运行 BLASTP
blastp -query "$pep2" -db "blastp/db" \
    -out "blastp/result.blast" \
    -evalue 1e-5 \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
    -num_threads $THREADS

cd blastp

cut -f1,2,4,3 result.blast > result.blast.res   # qid sid length pident

seqkit fx2tab -l -n -i "$cul/$pep1" > ref_protein_lengths.txt
seqkit fx2tab -l -n -i "$cul/$pep2" > query_protein_lengths.txt

python3 {SOFTWARE_HX}/Script/Homologous_gene_search/blast_merge.py \
    result.blast.res query_protein_lengths.txt ref_protein_lengths.txt blast.merge

# 筛选一致性差 或 覆盖率差的基因
awk -F'\t' 'BEGIN{OFS="\t"}
    $3=="gene" {
        match($9, /ID=([^;]+)/, a)
        if(a[1]!="") print a[1]
    }' "$cul/$gff1" > ../all_genes.list

awk -F "\t" '$3>50 && $7>0.5 && $8>0.5 {print $2}' blast.merge > ../blast_gene.tmp.list

cd "$cul"

grep -Fv -f blast_gene.tmp.list all_genes.list > all_new_genes.list

rm blast_gene.tmp.list

#=================== 提取基因所在位置 ===================
awk -F'\t' 'BEGIN{OFS="\t"}
    $3=="gene" {
        match($9, /ID=([^;]+)/, a)
        if(a[1]!="") print $1, $4-1, $5, a[1], $7
    }' "$gff1" > ref_genes.bed

awk -F'\t' 'BEGIN{OFS="\t"}
    $3=="gene" {
        match($9, /ID=([^;]+)/, a)
        if(a[1]!="") print $1, $4-1, $5, a[1], $7
    }' "$gff2" > query_genes.bed

grep -F -f all_new_genes.list ref_genes.bed > all_new_genes.bed

#=================== 分离 gap 区域基因 ===================
if [ -n "$gap" ] && [ -f "$gap" ]; then
    mkdir -p new_genes_Gap

    bedtools intersect \
        -a all_new_genes.bed \
        -b "$gap" \
        -wa \
        | awk '{print $4}' \
        | sort -u > new_genes_Gap/new_genes_Gap.list

    grep -F -f new_genes_Gap/new_genes_Gap.list ref_genes.bed \
        > new_genes_Gap/new_genes_Gap.bed

    grep -Fv -f new_genes_Gap/new_genes_Gap.list all_new_genes.list \
        > new_genes_none_Gap.list
else
    cp all_new_genes.list new_genes_none_Gap.list
fi

#=================== 获取非 gap 区的新注释基因序列 ===================
grep -F -f new_genes_none_Gap.list ref_genes.bed > new_genes_none_Gap.bed

bedtools getfasta \
    -fi "$genome1" \
    -bed new_genes_none_Gap.bed \
    -name -s \
    > new_genes_none_Gap.fa

#=================== BLASTN 到旧基因组 ===================
mkdir -p blastn

makeblastdb -in "$genome2" -dbtype nucl -out blastn/db

blastn -query new_genes_none_Gap.fa -db blastn/db \
    -out blastn/result.blast \
    -evalue 1e-5 \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
    -num_threads $THREADS

#=================== 区分其余新基因的类型 ===================
echo "[INFO] 开始筛选 blastn 最优匹配..."

mkdir -p best_hit
mkdir -p genes_anno_correct
mkdir -p new_genes_anno

# ------------------------------
# 1. 从 qseqid 中提取染色体编号（Chr1 / Chr2 ...）
# ------------------------------
extract_chr() {
    echo "$1" | sed -E 's/.*Chr([0-9]+).*/\1/'
}

# ------------------------------
# 2. 按优先级挑选每个 qseqid 的最佳匹配
# evalue 最小 > pident 最大 > length 最大
# 并要求染色体一致
# ------------------------------
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$9,$10,$11}' blastn/result.blast > blastn/tmp.rearranged


gawk -v OFS="\t" '
{
    q=$1; 
    s=$2; 
    pident=$3; 
    len=$4; 
    start_pos=$5; 
    end_pos=$6; 
    evalue=$7;

    # 提取 qseqid 染色体 Chr1
    chr_q = q; 
    sub(/.*Chr/, "", chr_q); 
    sub(/:.*/, "", chr_q);

    # 提取 sseqid 染色体 chr1
    chr_s = s; 
    sub(/.*chr/, "", chr_s);

    if (chr_q != chr_s) next;

    key = q;

    if (!(key in best)) {
        best[key]=$0;
        evalue_min[key]=evalue;
        pident_max[key]=pident;
        len_max[key]=len;
    } else {
        if (evalue < evalue_min[key] ||
            (evalue == evalue_min[key] && pident > pident_max[key]) ||
            (evalue == evalue_min[key] && pident == pident_max[key] && len > len_max[key])) {

            best[key]=$0;
            evalue_min[key]=evalue;
            pident_max[key]=pident;
            len_max[key]=len;
        }
    }
}
END {
    for (k in best) print best[k];
}
' blastn/tmp.rearranged > best_hit/best_hit.raw

echo "[INFO] 提取最佳匹配完毕 → best_hit/best_hit.raw"

# ------------------------------
# 3. 只保留 qseqid, sseqid, sstart, send
# ------------------------------
awk -v OFS="\t" '{print $1, $2, $5, $6}' best_hit/best_hit.raw > best_hit/best_hit.simple

echo "[INFO] 最佳匹配简化文件 → best_hit/best_hit.simple"

# ------------------------------
# 4. 与 query_genes.bed 比对是否有重叠
# ------------------------------
# 将 best_hit.simple 转为 BED
awk -v OFS="\t" '
{
    chrom=$2;

    # sstart=$3; send=$4;
    sstart=$3; 
    send=$4;

    # BED 要求 start < end
    if (sstart <= send) {
        start = sstart - 1;
        end = send;
    } else {
        start = send - 1;
        end = sstart;
    }

    print chrom, start, end, $1;
}
' best_hit/best_hit.simple > best_hit/best_hit.bed

# 进行 overlap 检查
bedtools intersect -wa -wb \
    -a best_hit/best_hit.bed \
    -b query_genes.bed \
    > best_hit/overlap.txt

# ------------------------------
# 5. 输出分类结果
# ------------------------------
# 获取修正基因
# 有重叠 → 输出 qseqid 与 query gene
awk -v OFS="\t" '
{
    qseq=$4;
    query_gene=$8;
    print qseq, query_gene;
}
' best_hit/overlap.txt \
> genes_anno_correct/genes_anno_correct.tsv

# 清洗基因名称
awk -v OFS="\t" '{
    gsub(/::.*/, "", $1);
    print $1, $2;
}' genes_anno_correct/genes_anno_correct.tsv > genes_anno_correct/genes_anno_correct.cleaned.tsv \
    && mv genes_anno_correct/genes_anno_correct.cleaned.tsv genes_anno_correct/genes_anno_correct.tsv

cut -f1 genes_anno_correct/genes_anno_correct.tsv | sort -u > genes_anno_correct/genes_anno_correct.list

grep -F -f genes_anno_correct/genes_anno_correct.list ref_genes.bed \
    > genes_anno_correct/genes_anno_correct.bed

echo "[INFO] 修正的基因 → genes_anno_correct/genes_anno_correct.list"

# 获取新注释的基因
# 无重叠 → 输出 qseqid 列表
grep -Fv -f genes_anno_correct/genes_anno_correct.list new_genes_none_Gap.list > new_genes_anno/new_genes_anno.list

grep -F -f new_genes_anno/new_genes_anno.list ref_genes.bed \
    > new_genes_anno/new_genes_anno.bed

echo "[INFO] 新注释的基因 → new_genes_anno/new_genes_anno.list"

echo "[INFO] 任务完成！"
