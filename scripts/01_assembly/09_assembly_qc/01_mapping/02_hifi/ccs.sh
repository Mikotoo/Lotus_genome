#!/bin/bash
#CSUB -J ccs
#CSUB -q c01
#CSUB -o ccs.out
#CSUB -e ccs.error
#CSUB -n 64
#CSUB -R span[hosts=1]

#minimap2 -ax map-pb {PROJ_LOTUS_ZC}/Gifu/00.assemble_result/Gifu_v1.0.fasta {PROJ_LOTUS_ZC}/data/Gifu/hifi/Gifu_hifi.fastq.gz -o ccs.sam
#samtools flagstat ccs.sam > ccs.flagstat
#samtools sort -@ 64 -m 10G -o ccs.bam -T ccs_tmp.ali ccs.sam
#samtools flagstat ccs.bam > ccs.flagstat
#samtools index ccs.bam

#samtools view -h -F 2304 ccs.bam -o primary_only.bam
#samtools view -b -q 30 -o ccs_q30.bam ccs.bam

#bedtools coverage -a {PROJ_04_LOTUS_GENOME}/00_assembly_result/GifuT2T/Gifu_100k.bed -b ccs_q30.bam -mean > gifu_hifi_100k_coverage.txt
#bedtools coverage -a {PROJ_04_LOTUS_GENOME}/00_assembly_result/GifuT2T/Gifu_10k.bed -b ccs_q30.bam -mean > gifu_hifi_10k_coverage.txt