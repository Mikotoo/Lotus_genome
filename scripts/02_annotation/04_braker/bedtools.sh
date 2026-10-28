bedtools intersect -a gene_filter.bed -b gene_filter.bed -wa -wb \
| awk -F'\t' '$4 < $9 {print $0}' \
> overlaps.tsv
