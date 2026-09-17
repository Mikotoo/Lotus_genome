#!/bin/bash
#CSUB -J process
#CSUB -q c01
#CSUB -o process.out
#CSUB -e process.error
#CSUB -n 88
#CSUB -R span[hosts=1]

bedtools getfasta -fi {HOME_ZC}/Project/Ref/ChiHei/ChiHei_v1.0.fasta -bed Chr13.bed > Chr13.fa
samtools faidx Chr13.fa
bedtools makewindows -g Chr13.fa.fai -w 2000 > Chr13_output.bed
bedtools getfasta -fi Chr13.fa -bed Chr13_output.bed > Chr13_output.2000.fasta

minimap2 -f 1000 -s 400 -ax ava-ont -d Chr13_output.fasta.mmi Chr13_output.2000.fasta
minimap2 -t 88 -f 10000 -s 400 -ax ava-ont --dual=yes --eqx Chr13_output.fasta.mmi Chr13_output.2000.fasta | samtools sort -m 4G -o Chr13_output.bam
python ~/software/StainedGlass-0.6/workflow/scripts/samIdentity.py --threads 88 --matches 400 --header Chr13_output.bam > Chr13_output.tbl
bgzip -c Chr13_output.tbl > Chr13_output.tbl.gz
python ~/software/StainedGlass-0.6/workflow/scripts/refmt.py --window 2000 --fai Chr13.fa.fai --full Chr13_output.full.tbl.gz Chr13_output.tbl.gz Chr13_output.bed.gz
mkdir -p results/Chr13_figures/pdfs
mkdir -p results/Chr13_figures/pngs
Rscript ~/software/StainedGlass-0.6/workflow/scripts/aln_plot.R -b Chr13_output.bed.gz --threads 88 --prefix Chr13



bedtools getfasta -fi {HOME_ZC}/Project/Ref/ChiHei/ChiHei_v1.0.fasta -bed Chr19.bed > Chr19.fa
samtools faidx Chr19.fa
bedtools makewindows -g Chr19.fa.fai -w 2000 > Chr19_output.bed
bedtools getfasta -fi Chr19.fa -bed Chr19_output.bed > Chr19_output.2000.fasta

minimap2 -f 1000 -s 400 -ax ava-ont -d Chr19_output.fasta.mmi Chr19_output.2000.fasta
minimap2 -t 88 -f 10000 -s 400 -ax ava-ont --dual=yes --eqx Chr19_output.fasta.mmi Chr19_output.2000.fasta | samtools sort -m 4G -o Chr19_output.bam
python ~/software/StainedGlass-0.6/workflow/scripts/samIdentity.py --threads 88 --matches 400 --header Chr19_output.bam > Chr19_output.tbl
bgzip -c Chr19_output.tbl > Chr19_output.tbl.gz
python ~/software/StainedGlass-0.6/workflow/scripts/refmt.py --window 2000 --fai Chr19.fa.fai --full Chr19_output.full.tbl.gz Chr19_output.tbl.gz Chr19_output.bed.gz
mkdir -p results/Chr19_figures/pdfs
mkdir -p results/Chr19_figures/pngs
Rscript ~/software/StainedGlass-0.6/workflow/scripts/aln_plot.R -b Chr19_output.bed.gz --threads 88 --prefix Chr19