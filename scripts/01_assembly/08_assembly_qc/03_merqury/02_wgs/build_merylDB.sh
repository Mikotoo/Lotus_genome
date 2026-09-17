#!/bin/bash
#CSUB -J meryl
#CSUB -q c01
#CSUB -o meryl.out
#CSUB -e meryl.error
#CSUB -n 88
#CSUB -R span[hosts=1]

meryl count k=21 {PROJ_LOTUS_ZC}/data/001-NC_raw_data/clean/Gifu_1-R1.fastq.gz output wgs1.merylDB
meryl count k=21 {PROJ_LOTUS_ZC}/data/001-NC_raw_data/clean/Gifu_1-R2.fastq.gz output wgs2.merylDB

meryl  union-sum output wgs.meryl wgs1.merylDB wgs2.merylDB
merqury.sh wgs.meryl {PROJ_LOTUS_ZC}/Gifu/00.assemble_result/Gifu_v1.0.fasta Gifu_wgsmer