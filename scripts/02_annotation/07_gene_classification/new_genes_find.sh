#!/bin/bash
#CSUB -J FindNewGenes
#CSUB -q c01
#CSUB -o FindNewGenes.out
#CSUB -e FindNewGenes.error
#CSUB -n 64
#CSUB -R span[hosts=1]

sh new_genes_finder.sh Lotus_GifuT2T.pep.fa\
    Lotus_GifuT2T_v1.0.fasta Lotus_GifuT2T.gff3\
    Lotus_Gifu_Old.pep.fa\
    Lotus_Gifu_Old.fa \
    Lotus_Gifu_Old.gff3\
    Old_region.uncovered.bed

sh new_genes_finder.sh Lotus_MG20T2T.pep.fa\
    Lotus_MG20T2T_v1.0.fasta Lotus_MG20T2T.gff3\
    Lotus_MG20_Old.pep.fa\
    Lotus_MG20_Old.fa \
    Lotus_MG20_Old.gff3\
    Old_region.uncovered.bed