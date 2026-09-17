#!/bin/bash
#CSUB -J fastp
#CSUB -q c01
#CSUB -o fastp.out
#CSUB -e fastp.error
#CSUB -n 88
#CSUB -R span[hosts=1]

cat align.conf|while read id;
do
arr=($id)
sample=${arr[0]}
fq1=${arr[1]}
fq2=${arr[2]}

fastp -i $fq1 -I $fq2 -o clean/${sample}-R1.fastq.gz -O clean/${sample}-R2.fastq.gz -w 88 -j report/${sample}.json -h report/${sample}.html

done