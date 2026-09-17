#!/bin/bash
#CSUB -J process
#CSUB -q c01
#CSUB -o process.out
#CSUB -e process.error
#CSUB -n 64
#CSUB -R span[hosts=1]

minimap2 -ax asm20 -t 64 --eqx {PROJ_LOTUS_ZC}/comp/reference/Gifu_v1.fa \
{PROJ_LOTUS_ZC}/comp/reference/Gifu_T2T.fa \
| samtools sort -O BAM - > Gifu_v1_Gifu_T2T.bam
samtools index Gifu_v1_Gifu_T2T.bam

syri -c Gifu_v1_Gifu_T2T.bam -r {PROJ_LOTUS_ZC}/comp/reference/Gifu_v1.fa \
-q {PROJ_LOTUS_ZC}/comp/reference/Gifu_T2T.fa -F B --prefix Gifu_v1_Gifu_T2T


minimap2 -ax asm20 -t 64 --eqx {PROJ_LOTUS_ZC}/comp/reference/Gifu_T2T.fa \
{PROJ_LOTUS_ZC}/comp/reference/MG20_T2T.fa \
| samtools sort -O BAM - > Gifu_T2T_MG20_T2T.bam
samtools index Gifu_T2T_MG20_T2T.bam

syri -c Gifu_T2T_MG20_T2T.bam -r {PROJ_LOTUS_ZC}/comp/reference/Gifu_T2T.fa \
-q {PROJ_LOTUS_ZC}/comp/reference/MG20_T2T.fa -F B --prefix Gifu_T2T_MG20_T2T


minimap2 -ax asm20 -t 64 --eqx {PROJ_LOTUS_ZC}/comp/reference/MG20_T2T.fa \
{PROJ_LOTUS_ZC}/comp/reference/MG20_v1.fa \
| samtools sort -O BAM - > MG20_T2T_MG20_v1.bam
samtools index MG20_T2T_MG20_v1.bam

syri -c MG20_T2T_MG20_v1.bam -r {PROJ_LOTUS_ZC}/comp/reference/MG20_T2T.fa \
-q {PROJ_LOTUS_ZC}/comp/reference/MG20_v1.fa -F B --prefix MG20_T2T_MG20_v1


plotsr --sr Gifu_v1_Gifu_T2Tsyri.out \
--sr Gifu_T2T_MG20_T2Tsyri.out \
--sr MG20_T2T_MG20_v1syri.out \
--genomes genomes.txt \
-o syri_samechromosome_plot.pdf