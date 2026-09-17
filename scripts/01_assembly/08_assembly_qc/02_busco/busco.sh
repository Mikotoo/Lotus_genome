#!/bin/bash
#CSUB -J busco
#CSUB -q c01
#CSUB -o busco.out
#CSUB -e busco.error
#CSUB -n 64
#CSUB -R span[hosts=1]

busco -i {PROJ_LOTUS_ZC}/Gifu/00.assemble_result/Gifu_v1.0.fasta \
      -l {SOFTWARE_ZC}/busco_downloads/lineages/embryophyta_odb10 \
      -o Gifu_lotus \
      -m genome \
      --offline --cpu 64