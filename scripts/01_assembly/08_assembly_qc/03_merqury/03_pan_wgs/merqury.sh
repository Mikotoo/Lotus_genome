#!/bin/bash
#CSUB -J merqury
#CSUB -q c01
#CSUB -o merqury.out
#CSUB -e merqury.error
#CSUB -n 64
#CSUB -R span[hosts=1]

cat {PROJ_LOTUS_ZC}/data/001-NC_raw_data/clean/align.conf|while read id;
do
arr=($id)
sample=${arr[0]}
fq1=${arr[1]}
fq2=${arr[2]}

cd {PROJ_LOTUS_ZC}/Gifu/08.quality/03.merqury/03.pan_wgsmer
mkdir $sample && cd $sample
meryl count k=21 $fq1 output ${sample}_R1.merylDB
meryl count k=21 $fq2 output ${sample}_R2.merylDB
meryl union-sum output ${sample}.meryl ${sample}_R1.merylDB ${sample}_R2.merylDB

merqury.sh ${sample}.meryl {PROJ_LOTUS_ZC}/Gifu/00.assemble_result/Gifu_v1.0.fasta Gifu_${sample}_wgsmer
done