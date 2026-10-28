#4.2 初步检查每条contig是否组装到端粒，辅助后续挂载排序

python get_tel.py
sort -u Gifu.tel_unit.bed |sort -k 1,1 -k 1n > tmp && mv tmp Gifu.tel_unit.bed
samtools faidx contigs.fa
cat contigs.fa.fai |awk '{print $1"\t"$2}' > Gifu.length
bedtools makewindows -g Gifu.length -w 1000 > Gifu_1k.bed
bedtools intersect -a Gifu_1k.bed -b Gifu.tel_unit.bed -wa |uniq -c|awk '$1>50{print $2"\t"$3"\t"$4"\t"$1}' > Gifu_tel.result