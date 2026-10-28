#!/bin/bash
#CSUB -J blast
#CSUB -q c01
#CSUB -o blast.out
#CSUB -e blast.error
#CSUB -n 88
#CSUB -R span[hosts=1]

set -euo pipefail

#=================== 参数设置 ===================
THREADS=88
cul=$(pwd)

#=================== 输入文件 ===================
GifuOld="Lotus_GifuOld.pep.fa"
GifuT2T="Lotus_GifuT2T.pep.fa"
MG20Old="Lotus_MG20Old.pep.fa"
MG20T2T="Lotus_MG20T2T.pep.fa"

ALL_GENE_LIST="{PROJ_LOTUS_HX}/03_annotation/Ortholog/04.blast_Ortholog3/GifuT2T_all_gene_list.txt"

#=================== 函数定义 ===================
run_blastp() {
    query=$1
    subject=$2
    outdir=$3

    mkdir -p "$outdir"

    makeblastdb -in "$subject" -dbtype prot -out "$outdir/db"

    blastp -query "$query" -db "$outdir/db" \
        -out "$outdir/result.blast" \
        -evalue 1e-5 \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
        -num_threads $THREADS
}

#=================== BLAST 比对 ===================
run_blastp "$GifuT2T" "$GifuOld"  "GifuT2T_GifuOld"
run_blastp "$GifuT2T" "$MG20Old"  "GifuT2T_MG20Old"
run_blastp "$GifuT2T" "$MG20T2T"  "GifuT2T_MG20T2T"

echo "✅ BLASTP 完成"

#=================== 结果处理 ===================
find ./GifuT2T_* -type f -name "*.blast" | while read blast; do

    blast1=$(basename "$blast")
    blast2=$(basename "$(dirname "$blast")")

    ref=${blast2##*_}
    query=${blast2%%_*}

    workdir="$cul/$blast2"
    cd "$workdir"

    # ------------------------------------------------
    # step0: 基础整理
    # ------------------------------------------------
    cut -f1,2,4,3 "$blast1" > "${blast1}.res"

    seqkit fx2tab -l -n -i "$cul/Lotus_${ref}.pep.fa"   > ref_len.txt
    seqkit fx2tab -l -n -i "$cul/Lotus_${query}.pep.fa" > query_len.txt

    python3 $cul/blast_merge.py \
        "${blast1}.res" \
        query_len.txt \
        ref_len.txt \
        blast1.merge

    # ------------------------------------------------
    # 1️⃣ 去版本号 + 基础过滤 → tmp1
    # ------------------------------------------------
    awk -F'\t' '
    BEGIN{OFS="\t"}
    {
        gsub(/\.[0-9]+$/, "", $1)
        gsub(/\.[0-9]+$/, "", $2)
        if ($3 > 50)
            print
    }' blast1.merge > tmp1.tsv

    # ------------------------------------------------
    # 2️⃣ 提取 sure pairs
    #   100% identity + 完全覆盖
    # ------------------------------------------------
    awk -F'\t' '
    BEGIN{OFS="\t"}
    $3==100 && $7==1 && $8==1 {print $1,$2}
    ' tmp1.tsv | sort -u > sure_pairs.tsv

    # ------------------------------------------------
    # 3️⃣ 去除已是 sure 的 query → tmp2
    # ------------------------------------------------
    awk -F'\t' '
    NR==FNR {a[$1]; next}
    !($1 in a)
    ' sure_pairs.tsv tmp1.tsv > tmp2.tsv

    # ------------------------------------------------
    # 4️⃣ 对剩余 query，按 $9 最大选 best hit
    # ------------------------------------------------
    awk -F'\t' '
    BEGIN{OFS="\t"}
    {
        if (!($1 in max) || $9 > max[$1]) {
            max[$1] = $9
            line[$1] = $1 OFS $2
        }
    }
    END{
        for (i in line) print line[i]
    }' tmp2.tsv > best_pairs.tsv

    # ------------------------------------------------
    # 5️⃣ 合并 sure + best
    # ------------------------------------------------
    cat sure_pairs.tsv best_pairs.tsv | sort -u > all_pairs.tsv

    # ------------------------------------------------
    # 6️⃣ 补全所有 GifuT2T 基因
    # ------------------------------------------------
    awk -F'\t' '
    BEGIN{OFS="\t"}
    NR==FNR {map[$1]=$2; next}
    {
        if ($1 in map)
            print $1, map[$1]
        else
            print $1, "NA"
    }' all_pairs.tsv "$ALL_GENE_LIST" > final_pairs.tsv

    echo "✅ [$blast2] finished → final_pairs.tsv"

    cd "$cul"
done

echo "🎉 所有比对与基因配对流程完成"

#================================================================
# 7️⃣ 合并所有 final_pairs.tsv → 四列表
#================================================================
echo "🔗 合并所有 final_pairs.tsv → All_species_merged_by_GifuT2T.tsv"

paste \
    "$cul/GifuT2T_GifuOld/final_pairs.tsv" \
    "$cul/GifuT2T_MG20Old/final_pairs.tsv" \
    "$cul/GifuT2T_MG20T2T/final_pairs.tsv" \
| awk -F'\t' '
BEGIN{
    OFS="\t";
    print "GifuT2T","GifuOld","MG20Old","MG20T2T"
}
{
    # 每个 final_pairs.tsv 都是两列：GifuT2T  Other
    # paste 后为：1 2 | 3 4 | 5 6
    print $1, $2, $4, $6
}
' > All_species_merged_by_GifuT2T.tsv

echo "🎉 所有物种合并完成 → All_species_merged_by_GifuT2T.tsv"