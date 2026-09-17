#!/bin/bash
#CSUB -J bwa
#CSUB -q c01
#CSUB -o bwa.out
#CSUB -e bwa.error
#CSUB -n 88
#CSUB -R span[hosts=1]

#bwa index -p Gifu_v1.0.fasta -a bwtsw Gifu_v1.0.fasta
#bwa mem Gifu_v1.0.fasta {PROJ_LOTUS_ZC}/data/001-NC_raw_data/clean/Gifu_1-R1.fastq.gz {PROJ_LOTUS_ZC}/data/001-NC_raw_data/clean/Gifu_1-R2.fastq.gz > bwa_1.sam
#samtools flagstat bwa_1.sam > bwa_1.flagstat

#bwa mem Gifu_v1.0.fasta {PROJ_LOTUS_ZC}/data/001-NC_raw_data/clean/Gifu_2-R1.fastq.gz {PROJ_LOTUS_ZC}/data/001-NC_raw_data/clean/Gifu_2-R2.fastq.gz > bwa_2.sam
#samtools flagstat bwa_2.sam > bwa_2.flagstat
#samtools sort -@ 88 -o bwa.bam -T bwa_tmp.ali bwa.sam

samtools view -b -q 30 -o wgs_q30.bam {PROJ_04_LOTUS_GENOME}/06_SNP/02.toCram/picard/GifuT2T_markdup.bam


bedtools coverage -a {PROJ_04_LOTUS_GENOME}/00_assembly_result/GifuT2T/Gifu_100k.bed -b wgs_q30.bam -mean > gifu_wgs_100k_coverage.txt
bedtools coverage -a {PROJ_04_LOTUS_GENOME}/00_assembly_result/GifuT2T/Gifu_10k.bed -b wgs_q30.bam -mean > gifu_wgs_10k_coverage.txt