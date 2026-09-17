bedtools intersect -a gene.bed -b gene.bed -wa -wb \
| awk -F'\t' '$4 < $8 {print $0}' \
> overlaps.tsv
