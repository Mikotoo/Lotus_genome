#!/bin/bash
#CSUB -J merqury
#CSUB -q c01
#CSUB -o merqury.out
#CSUB -e merqury.error
#CSUB -n 64
#CSUB -R span[hosts=1]

merqury.sh wgs.meryl {PROJ_LOTUS_ZC}/Gifu/00.assemble_result/Gifu_v1.0.fasta Gifu_wgsmer