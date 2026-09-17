bedtools intersect -a Lotus_Gifu_gene.bed -b Lotus_Gifu_gene.bed -wa -wb \
| awk -F'\t' '$4 < $9 {print $0}' \
> overlaps.tsv
