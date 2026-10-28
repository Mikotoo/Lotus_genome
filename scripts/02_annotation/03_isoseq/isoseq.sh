#!/bin/bash
#CSUB -J isoseq
#CSUB -q c02
#CSUB -o isoseq.out
#CSUB -e isoseq.error
#CSUB -n 88
#CSUB -R span[ptile=88]

sample1={PROJ_LOTUS_ZC}/data/Gifu/isoseq/hifi_reads/P22TR251406648-1-bc12_r84069_20250617_081306_1_A01.hifi_reads.bam
primer1={PROJ_LOTUS_ZC}/data/Gifu/isoseq/hifi_reads/P22TR251406648-1-bc12_r84069_20250617_081306_1_A01.primer.fasta

isoseq3 refine "$sample1" "$primer1" Gifu.refined.bam
isoseq3 cluster Gifu.refined.bam Gifu.clustered.bam --verbose --use-qvs


samtools fastq Gifu.clustered.bam  > Gifu.clustered.fastq
minimap2 -ax splice -uf -k14 {PROJ_LOTUS_ZC}/Gifu/00.assemble_result/Gifu_v1.0.fasta Gifu.clustered.fastq > Gifu_isoforms.sam
samtools sort -o Gifu_isoforms.sorted.bam Gifu_isoforms.sam
samtools index Gifu_isoforms.sorted.bam
