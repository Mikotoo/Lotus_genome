#!/bin/bash
#CSUB -J filter
#CSUB -q c01
#CSUB -o filter.out
#CSUB -e filter.error
#CSUB -n 88
#CSUB -R span[hosts=1]

python AnnotationFilter.py -l 20 \
-i braker.aa -o Gifu_filtered.pep.fa \
-c braker.codingseq -C Gifu_filtered.cds.fa \
-g braker.gff3 -G Gifu_filtered.gff3