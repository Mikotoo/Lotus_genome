#!/bin/bash
#CSUB -J process
#CSUB -q c01
#CSUB -o process.out
#CSUB -e process.error
#CSUB -n 88
#CSUB -R span[hosts=1]

#python3 canon_unit_and_length.py \
#  -i mg20_cent.bed \
#  -o mg20.units.canonical.tsv \
#  --summary mg20.length_summary.tsv \
#  --seq-col 9

#####
#in="gifu.units.canonical.tsv"   # 改成你的文件名
#out="mg20_chr_length_summary.tsv"

#awk -F'\t' 'NR>1{
#  chr=$1; L=$5; span=$8;
#  key=chr"\t"L;
#  cov[key]+=span; cnt[key]+=1
#}
#END{
#  print "chr\tunit_len\tcount\ttotal_span_bp";
#  for(k in cov) print k"\t"cnt[k]"\t"cov[k]
#}' "$in" | sort -k1,1 -k4,4nr > "$out"


#每条染色体覆盖最大的长度

#awk -F'\t' 'NR>1{key=$1"\t"$5; cov[key]+=$8}
#END{for(k in cov) print k"\t"cov[k]}' "$in" \
#| sort -k1,1 -k3,3nr \
#| awk 'BEGIN{FS=OFS="\t"} {if($1!=p){print; p=$1}}'

#in="mg20.units.canonical.tsv"
#out="mg20_chr_top_monomer.tsv"

#awk -F'\t' 'NR>1{
#  chr=$1; L=$5; unit=$6; span=$8;
#  if(L>=100 && L<=400){
#    key=chr"\t"unit;
#    cov[key]+=span;
#    len[key]=L;
#  }
#}
#END{
#  for(k in cov){
#    split(k,a,"\t");
#    chr=a[1]; unit=a[2];
#    print chr"\t"len[k]"\t"unit"\t"cov[k];
#  }
#}' "$in" | sort -k1,1 -k4,4nr \
#| awk 'BEGIN{FS=OFS="\t"}
#{
#  chr=$1;
#  top[chr]++
#  if(top[chr]<=10) print
#}' > "$out"

cat mg20_chr_top_monomer.tsv gifu_chr_top_monomer.tsv > lotus_chr_top_monomer.tsv

in="lotus_chr_top_monomer.tsv"         # 你的这个输出文件名
fa="lotus_chr_top_monomer.fa"

awk 'BEGIN{FS=OFS="\t"}
{
  chr=$1; L=$2; seq=$3; cov=$4;
  id=chr"|L="L"|cov="cov"|n="NR;
  print ">"id"\n"seq
}' "$in" > "$fa"


makeblastdb -in "$fa" -dbtype nucl -out lotus_monomer
blastn -task blastn-short \
  -query "$fa" -db lotus_monomer \
  -dust no -soft_masking false \
  -word_size 7 \
  -evalue 1e-10 \
  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
  -max_target_seqs 10000 \
  > self.allvsall.tsv

awk -F'\t' '$1!=$2' self.allvsall.tsv > self.noSelf.tsv
awk -F'\t' '($3>=80 && $4>=100)' self.noSelf.tsv > self.sim80.tsv
awk -F'\t' '($3>=90 && $4>=140)' self.noSelf.tsv > self.sim90.tsv
sort -k1,1 -k12,12nr self.sim90.tsv \
| awk -F'\t' '!seen[$1]++{print $1"\t"$2"\tpid="$3"\taln="$4"\tbits="$12}' \
> top_hit_per_query.tsv

awk -F'\t' '{
  split($1,a,"|"); split($2,b,"|");
  chr1=a[1]; chr2=b[1];
  key=chr1"\t"chr2;
  n[key]++
} END{
  for(k in n) print k"\t"n[k]
}' self.sim90.tsv | sort -k3,3nr > chr_pair_links.tsv


