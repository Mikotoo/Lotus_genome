#!/bin/bash
#CSUB -J meryl
#CSUB -q c01
#CSUB -o meryl.out
#CSUB -e meryl.error
#CSUB -n 88
#CSUB -R span[hosts=1]

meryl count k=21 {PROJ_LOTUS_ZC}/data/Gifu/hifi/Gifu_hifi.fastq.gz output ccs.merylDB

