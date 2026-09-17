#!/bin/bash
#CSUB -J hifiasm
#CSUB -q c01
#CSUB -o hifiasm.out
#CSUB -e hifiasm.error
#CSUB -n 88
#CSUB -R span[hosts=1]

#hifiasm -t 88 -o Gifu.sam --telo-m CCCTAAA --ul {PROJ_LOTUS_ZC}/data/Gifu/ont/Gifu-jing-ye.pass.ul.fq.gz {PROJ_LOTUS_ZC}/data/Gifu/hifi/Gifu_hifi.fastq.gz

gfatools gfa2fa Gifu.sam.bp.p_ctg.gfa > Gifu.hifiasm.fasta
samtools faidx Gifu.hifiasm.fasta