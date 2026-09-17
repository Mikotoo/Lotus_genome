#!/bin/bash
#CSUB -J process
#CSUB -q c02
#CSUB -o process.out
#CSUB -e process.error
#CSUB -n 64
#CSUB -R span[hosts=1]

#makeblastdb -in Lotus_GifuT2T_v1.0.fasta -dbtype nucl -out gifu_genome
#makeblastdb -in Lotus_MG20T2T_v1.0.fasta -dbtype nucl -out mg20_genome

#blastn -task blastn-short \
#  -query cent_unit.fa -db gifu_genome \
#  -dust no -soft_masking false \
#  -word_size 7 \
#  -evalue 1e-10 \
#  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
#  > gifu_cent.blast

#blastn -task blastn-short \
#  -query cent_unit.fa -db mg20_genome \
#  -dust no -soft_masking false \
#  -word_size 7 \
#  -evalue 1e-10 \
#  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
#  > mg20_cent.blast


#blast=gifu_cent.blast

# 阈值可按需要调：例如 identity>=90 且 aln_len>=120
#PID=80
#ALN=100

#awk -v PID="$PID" -v ALN="$ALN" 'BEGIN{FS=OFS="\t"}
#{
#  q=$1; chr=$2; pid=$3; alen=$4;
#  s1=$9; s2=$10;

#  if(pid < PID) next;
#  if(alen < ALN) next;

  # strand: sstart < send -> + ; else -
#  strand = (s1 < s2 ? "+" : "-");

  # BED coordinates (0-based, half-open)
#  st = (s1 < s2 ? s1 : s2);
#  en = (s1 > s2 ? s1 : s2);

  # output BED6: chr start end name score strand
#  print chr, st-1, en, q, pid, strand
#}' "$blast" \
#| sort -k1,1 -k2,2n > gifu_units.all.bed

# 按 unit 拆成 4 个 bed
#for u in unit_161 unit_172 unit_186 unit_330; do
#  awk -v u=$u '$4==u' gifu_units.all.bed > gifu_${u}.bed
#done


#for a in mg20_unit_161 mg20_unit_172 mg20_unit_186 mg20_unit_330; do
#  for b in mg20_unit_161 mg20_unit_172 mg20_unit_186 mg20_unit_330; do
#    if [[ "$a" < "$b" ]]; then
#      echo "== $a vs $b (ROI) =="
#      bedtools jaccard -a ${a}.bed -b ${b}.bed
#    fi
#  done
#done


# 定义窗口
#echo -e "Chr5\t34000000\t34100000" > chr5_34.0_34.1.bed

# 提取窗口内命中
#for u in gifu_unit_161 gifu_unit_172 gifu_unit_186 gifu_unit_330; do
#  bedtools intersect -a ${u}.bed -b chr5_34.0_34.1.bed > ${u}.Chr5_34_34.1.bed
#done

# 可选：合并碎片，减少视觉“交错假象”
#for u in gifu_unit_161 gifu_unit_172 gifu_unit_186 gifu_unit_330; do
#  bedtools sort -i ${u}.Chr5_34_34.1.bed \
#  | bedtools merge -i - > ${u}.Chr5_34_34.1.merged.bed
#done


bash run_unit_mafft_consensus.sh
