#!/bin/bash
#CSUB -J purge
#CSUB -q c01
#CSUB -o purge.out
#CSUB -e purge.error
#CSUB -n 88
#CSUB -R span[hosts=1]

minimap2 -ax map-pb {PROJ_LOTUS_ZC}/Gifu/01.hifiasm/Gifu.hifiasm.fasta {PROJ_LOTUS_ZC}/data/Gifu/hifi/Gifu_hifi.fastq.gz \
|samtools view -hF 256 - \
|samtools sort -@ 20 -m 1G -o aligned.bam -T tmp.ali
