#!/bin/bash
#CSUB -J findmarker
#CSUB -q c01
#CSUB -o findmarker.out
#CSUB -e findmarker.error
#CSUB -n 64
#CSUB -R span[hosts=1]

# 输入文件
csv=./04_Anno/cluster_06_cluster_markers_gene.csv
fasta={PROJ_LOTUS_HX}/00_assembly_result/GifuT2T/Lotus_GifuT2T.pep.fa
annot=./Lotus_Gifu.emapper.annotations
outdir=./04_Anno

mkdir -p $outdir

# 需要处理的 cluster
#clusters=("3" "4" "6")
clusters=("8")

for n in "${clusters[@]}"; do
    echo ">>> Processing cluster $n ..."

    # Step1: 获取 top10 基因列表
    awk -F',' -v c="$n" '
        NR>1 {
            gsub(/"/, "", $0);                  
            if ($6 == c) {
                gsub(/-/, "_", $7);             
                print $2 "\t" $7                
            }
        }' "$csv" \
        | sort -k1,1nr \
        | head -n 20 \
        | cut -f2 \
        > ${outdir}/cluster_${n}_DEGs.list

    # Step2: 从 pep.fa 提取序列
    seqkit grep -f ${outdir}/cluster_${n}_DEGs.list "$fasta" \
        > ${outdir}/cluster_${n}_DEGs.pep.fa

    echo "Sequence extraction for cluster $n done."

    # Step3: 从 EggNOG 注释中提取对应基因的注释
    # 注释文件的第一列是 query ID
    grep -v "^#" "$annot" | awk -v OFS="\t" '
        NR==FNR {a[$1]=1; next}
        ($1 in a) {print}
    ' ${outdir}/cluster_${n}_DEGs.list - \
        > ${outdir}/cluster_${n}_DEGs.annotation

    echo "Annotation extraction for cluster $n done."
done

echo "All clusters finished."
