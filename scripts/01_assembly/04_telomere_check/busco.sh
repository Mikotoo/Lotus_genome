#!/bin/bash
#CSUB -J busco
#CSUB -q c01
#CSUB -o busco.out
#CSUB -e busco.error
#CSUB -n 64
#CSUB -R span[hosts=1]

busco -i contigs.fa \
      -l {SOFTWARE_ZC}/busco_downloads/lineages/embryophyta_odb10 \
      -o Gifu_lotus \
      -m genome \
      --cpu 64