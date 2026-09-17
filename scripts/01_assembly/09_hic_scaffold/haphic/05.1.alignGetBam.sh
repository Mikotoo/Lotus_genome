#!/bin/bash
#CSUB -J align
#CSUB -q c01
#CSUB -o align.out
#CSUB -e align.error
#CSUB -n 88
#CSUB -R span[hosts=1]

# (1) Align Hi-C data to the assembly, remove PCR duplicates and filter out secondary and supplementary alignments
#bwa index asm.fa
#bwa mem -5SP -t 48 asm.fa {PROJ_LOTUS_ZC}/data/Gifu/hic/Gifu_hic_1.fastq {PROJ_LOTUS_ZC}/data/Gifu/hic/Gifu_hic_2.fastq | samblaster | samtools view - -@ 48 -S -h -b -F 3340 -o HiC.bam

# (2) Filter the alignments with MAPQ 1 (mapping quality ≥ 1) and NM 3 (edit distance < 3)
{SOFTWARE_ZC}/HapHiC-main/utils/filter_bam HiC.bam 1 --nm 3 --threads 48 | samtools view - -b -@ 48 -o HiC.filtered.bam