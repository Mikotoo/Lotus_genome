#!/bin/bash
#CSUB -J bwa
#CSUB -q c01
#CSUB -o index.out
#CSUB -e index.error
#CSUB -n 88
#CSUB -R span[hosts=1]

bowtie2-build -f Lotus_GifuT2T_v1.0.fasta --thread 88 Gifu_T2T
