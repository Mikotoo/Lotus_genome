#!/bin/bash
#CSUB -J findmarker
#CSUB -q c01
#CSUB -o findmarker.out
#CSUB -e findmarker.error
#CSUB -n 64
#CSUB -R span[hosts=1]

cul=`pwd`

# 输入文件
csv=./04_Anno/cluster_06_cluster_markers_gene.csv
fasta={PROJ_LOTUS_HX}/00_assembly_result/GifuT2T/Lotus_GifuT2T.pep.fa
out=./04_Anno/cluster_DEGs.pep.fa

# 对DEGs进行过滤，提取 gene 列，去掉表头并去重，并将-替换为_
awk -F',' 'NR>1 {
    gsub(/"/, "", $0);
    if ( $5 < 0.05) {
        gsub(/-/, "_", $7);
        print $7
    }
}' "$csv" | sort -u > ./04_Anno/cluster_DEGs.list

# 用 seqkit 从 pep.fa 中提取对应序列
seqkit grep -f ./04_Anno/cluster_DEGs.list "$fasta" > "$out"

echo "✅ Done! Extracted sequences saved to: $out"

outdir=./04_Anno/blast_Gmax
mkdir -p "$outdir"

# 构建blast数据库
makeblastdb -in "./04_Anno/Gmax_marker_genes.pep.fa" -dbtype prot -out "$outdir/db"

# 运行blastp
blastp -query "$out" -db "$outdir/db" \
    -out "$outdir/result.blast" \
    -evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
    -num_threads 32

cut -f1,2,4,3 $outdir/result.blast > $outdir/result.blast.res

seqkit fx2tab -l -n -i ./04_Anno/Gmax_marker_genes.pep.fa > $outdir/ref_protein_lengths.txt

seqkit fx2tab -l -n -i $out > $outdir/query_protein_lengths.txt

cd $outdir

python3 $cul/blast_merge.py result.blast.res query_protein_lengths.txt ref_protein_lengths.txt blast1.merge

awk '$3>50 && $7>0.5 && $8>0.5' blast1.merge | \
awk '{
    diff = (1-$7>=0?1-$7:$7-1) + (1-$8>=0?1-$8:$8-1)
    key = $2
    if (!(key in best) || $3 > best_id[key] || ($3 == best_id[key] && diff < best_diff[key])) {
        best[key] = $0
        best_id[key] = $3
        best_diff[key] = diff
    }
}
END {
    for (k in best) print best[k]
}' > blast1.filtered.best

cd $cul
outdir=./04_Anno/blast_Gm2.1
mkdir -p "$outdir"

# 构建blast数据库
makeblastdb -in "./04_Anno/Gm2.1_marker_genes.pep.fa" -dbtype prot -out "$outdir/db"

# 运行blastp
blastp -query "$out" -db "$outdir/db" \
    -out "$outdir/result.blast" \
    -evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
    -num_threads 32

cut -f1,2,4,3 $outdir/result.blast > $outdir/result.blast.res

seqkit fx2tab -l -n -i ./04_Anno/Gm2.1_marker_genes.pep.fa > $outdir/ref_protein_lengths.txt

seqkit fx2tab -l -n -i $out > $outdir/query_protein_lengths.txt

cd $outdir

python3 $cul/blast_merge.py result.blast.res query_protein_lengths.txt ref_protein_lengths.txt blast1.merge

awk '$3>50 && $7>0.5 && $8>0.5' blast1.merge | \
awk '{
    diff = (1-$7>=0?1-$7:$7-1) + (1-$8>=0?1-$8:$8-1)
    key = $2
    if (!(key in best) || $3 > best_id[key] || ($3 == best_id[key] && diff < best_diff[key])) {
        best[key] = $0
        best_id[key] = $3
        best_diff[key] = diff
    }
}
END {
    for (k in best) print best[k]
}' > blast1.filtered.best

cd $cul
outdir=./04_Anno/blast_GmZH
mkdir -p "$outdir"

# 构建blast数据库
makeblastdb -in "./04_Anno/GmZH_marker_genes.pep.fa" -dbtype prot -out "$outdir/db"

# 运行blastp
blastp -query "$out" -db "$outdir/db" \
    -out "$outdir/result.blast" \
    -evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
    -num_threads 32

cut -f1,2,4,3 $outdir/result.blast > $outdir/result.blast.res

seqkit fx2tab -l -n -i ./04_Anno/GmZH_marker_genes.pep.fa > $outdir/ref_protein_lengths.txt

seqkit fx2tab -l -n -i $out > $outdir/query_protein_lengths.txt

cd $outdir

python3 $cul/blast_merge.py result.blast.res query_protein_lengths.txt ref_protein_lengths.txt blast1.merge

awk '$3>50 && $7>0.5 && $8>0.5' blast1.merge | \
awk '{
    diff = (1-$7>=0?1-$7:$7-1) + (1-$8>=0?1-$8:$8-1)
    key = $2
    if (!(key in best) || $3 > best_id[key] || ($3 == best_id[key] && diff < best_diff[key])) {
        best[key] = $0
        best_id[key] = $3
        best_diff[key] = diff
    }
}
END {
    for (k in best) print best[k]
}' > blast1.filtered.best

cd $cul
outdir=./04_Anno/blast_At
mkdir -p "$outdir"

# 构建blast数据库
makeblastdb -in "./04_Anno/At_marker_genes.pep.fa" -dbtype prot -out "$outdir/db"

# 运行blastp
blastp -query "$out" -db "$outdir/db" \
    -out "$outdir/result.blast" \
    -evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
    -num_threads 32

cut -f1,2,4,3 $outdir/result.blast > $outdir/result.blast.res

seqkit fx2tab -l -n -i ./04_Anno/At_marker_genes.pep.fa > $outdir/ref_protein_lengths.txt

seqkit fx2tab -l -n -i $out > $outdir/query_protein_lengths.txt

cd $outdir

python3 $cul/blast_merge.py result.blast.res query_protein_lengths.txt ref_protein_lengths.txt blast1.merge

awk '$3>50 && $7>0.5 && $8>0.5' blast1.merge | \
awk '{
    diff = (1-$7>=0?1-$7:$7-1) + (1-$8>=0?1-$8:$8-1)
    key = $2
    if (!(key in best) || $3 > best_id[key] || ($3 == best_id[key] && diff < best_diff[key])) {
        best[key] = $0
        best_id[key] = $3
        best_diff[key] = diff
    }
}
END {
    for (k in best) print best[k]
}' > blast1.filtered.best