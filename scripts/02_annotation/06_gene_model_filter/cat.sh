#!/bin/bash
#CSUB -J cat
#CSUB -q c01
#CSUB -o cat.out
#CSUB -e cat.error
#CSUB -n 64
#CSUB -R span[host=1]

cat {PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/2|while read id;do cat ${id}>> Gifu_RNA_all_R2.fastq.gz; done
cat {PROJ_LOTUS_ZC}/data/Gifu/RNA/clean/1|while read id;do cat ${id}>> Gifu_RNA_all_R1.fastq.gz; done