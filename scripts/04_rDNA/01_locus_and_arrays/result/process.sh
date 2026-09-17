#!/bin/bash

bedtools getfasta -fi ../Lotus_GifuT2T_v1.0.fasta -bed Chr2_0_3Mbp.bed > Chr2_0_3Mbp.fa
samtools faidx Chr2_0_3Mbp.fa
bedtools makewindows -g Chr2_0_3Mbp.fa.fai -w 2000 > Chr2_0_3Mbp_output.bed
bedtools getfasta -fi Chr2_0_3Mbp.fa -bed Chr2_0_3Mbp_output.bed > Chr2_0_3Mbp_output.2000.fasta

minimap2 -f 1000 -s 400 -ax ava-ont -d Chr2_0_3Mbp_output.fasta.mmi Chr2_0_3Mbp_output.2000.fasta
minimap2 -t 88 -f 10000 -s 400 -ax ava-ont --dual=yes --eqx Chr2_0_3Mbp_output.fasta.mmi Chr2_0_3Mbp_output.2000.fasta | samtools sort -m 4G -o Chr2_0_3Mbp_output.bam
python3.13 {SOFTWARE_LOCAL}/StainedGlass-0.6/workflow/scripts/samIdentity.py --threads 88 --matches 400 --header Chr2_0_3Mbp_output.bam > Chr2_0_3Mbp_output.tbl
bgzip -c Chr2_0_3Mbp_output.tbl > Chr2_0_3Mbp_output.tbl.gz
python3.13 {SOFTWARE_LOCAL}/StainedGlass-0.6/workflow/scripts/refmt.py --window 2000 --fai Chr2_0_3Mbp.fa.fai --full Chr2_0_3Mbp_output.full.tbl.gz Chr2_0_3Mbp_output.tbl.gz Chr2_0_3Mbp_output.bed.gz
mkdir -p results/Chr2_0_3Mbp_figures/pdfs
mkdir -p results/Chr2_0_3Mbp_figures/pngs
Rscript {SOFTWARE_LOCAL}/StainedGlass-0.6/workflow/scripts/aln_plot.R -b Chr2_0_3Mbp_output.bed.gz --threads 88 --prefix Chr2_0_3Mbp



bedtools getfasta -fi ../Lotus_GifuT2T_v1.0.fasta -bed Chr2_24_26Mbp.bed > Chr2_24_26Mbp.fa
samtools faidx Chr2_24_26Mbp.fa
bedtools makewindows -g Chr2_24_26Mbp.fa.fai -w 2000 > Chr2_24_26Mbp_output.bed
bedtools getfasta -fi Chr2_24_26Mbp.fa -bed Chr2_24_26Mbp_output.bed > Chr2_24_26Mbp_output.2000.fasta

minimap2 -f 1000 -s 400 -ax ava-ont -d Chr2_24_26Mbp_output.fasta.mmi Chr2_24_26Mbp_output.2000.fasta
minimap2 -t 88 -f 10000 -s 400 -ax ava-ont --dual=yes --eqx Chr2_24_26Mbp_output.fasta.mmi Chr2_24_26Mbp_output.2000.fasta | samtools sort -m 4G -o Chr2_24_26Mbp_output.bam
python3.13 {SOFTWARE_LOCAL}/StainedGlass-0.6/workflow/scripts/samIdentity.py --threads 88 --matches 400 --header Chr2_24_26Mbp_output.bam > Chr2_24_26Mbp_output.tbl
bgzip -c Chr2_24_26Mbp_output.tbl > Chr2_24_26Mbp_output.tbl.gz
python3.13 {SOFTWARE_LOCAL}/StainedGlass-0.6/workflow/scripts/refmt.py --window 2000 --fai Chr2_24_26Mbp.fa.fai --full Chr2_24_26Mbp_output.full.tbl.gz Chr2_24_26Mbp_output.tbl.gz Chr2_24_26Mbp_output.bed.gz
mkdir -p results/Chr2_24_26Mbp_figures/pdfs
mkdir -p results/Chr2_24_26Mbp_figures/pngs
Rscript {SOFTWARE_LOCAL}/StainedGlass-0.6/workflow/scripts/aln_plot.R -b Chr2_24_26Mbp_output.bed.gz --threads 88 --prefix Chr2_24_26Mbp



bedtools getfasta -fi ../Lotus_GifuT2T_v1.0.fasta -bed Chr5_14_16Mbp.bed > Chr5_14_16Mbp.fa
samtools faidx Chr5_14_16Mbp.fa
bedtools makewindows -g Chr5_14_16Mbp.fa.fai -w 2000 > Chr5_14_16Mbp_output.bed
bedtools getfasta -fi Chr5_14_16Mbp.fa -bed Chr5_14_16Mbp_output.bed > Chr5_14_16Mbp_output.2000.fasta

minimap2 -f 1000 -s 400 -ax ava-ont -d Chr5_14_16Mbp_output.fasta.mmi Chr5_14_16Mbp_output.2000.fasta
minimap2 -t 88 -f 10000 -s 400 -ax ava-ont --dual=yes --eqx Chr5_14_16Mbp_output.fasta.mmi Chr5_14_16Mbp_output.2000.fasta | samtools sort -m 4G -o Chr5_14_16Mbp_output.bam
python3.13 {SOFTWARE_LOCAL}/StainedGlass-0.6/workflow/scripts/samIdentity.py --threads 88 --matches 400 --header Chr5_14_16Mbp_output.bam > Chr5_14_16Mbp_output.tbl
bgzip -c Chr5_14_16Mbp_output.tbl > Chr5_14_16Mbp_output.tbl.gz
python3.13 {SOFTWARE_LOCAL}/StainedGlass-0.6/workflow/scripts/refmt.py --window 2000 --fai Chr5_14_16Mbp.fa.fai --full Chr5_14_16Mbp_output.full.tbl.gz Chr5_14_16Mbp_output.tbl.gz Chr5_14_16Mbp_output.bed.gz
mkdir -p results/Chr5_14_16Mbp_figures/pdfs
mkdir -p results/Chr5_14_16Mbp_figures/pngs
Rscript {SOFTWARE_LOCAL}/StainedGlass-0.6/workflow/scripts/aln_plot.R -b Chr5_14_16Mbp_output.bed.gz --threads 88 --prefix Chr5_14_16Mbp



bedtools getfasta -fi ../Lotus_GifuT2T_v1.0.fasta -bed Chr6_8_10Mbp.bed > Chr6_8_10Mbp.fa
samtools faidx Chr6_8_10Mbp.fa
bedtools makewindows -g Chr6_8_10Mbp.fa.fai -w 2000 > Chr6_8_10Mbp_output.bed
bedtools getfasta -fi Chr6_8_10Mbp.fa -bed Chr6_8_10Mbp_output.bed > Chr6_8_10Mbp_output.2000.fasta

minimap2 -f 1000 -s 400 -ax ava-ont -d Chr6_8_10Mbp_output.fasta.mmi Chr6_8_10Mbp_output.2000.fasta
minimap2 -t 88 -f 10000 -s 400 -ax ava-ont --dual=yes --eqx Chr6_8_10Mbp_output.fasta.mmi Chr6_8_10Mbp_output.2000.fasta | samtools sort -m 4G -o Chr6_8_10Mbp_output.bam
python3.13 {SOFTWARE_LOCAL}/StainedGlass-0.6/workflow/scripts/samIdentity.py --threads 88 --matches 400 --header Chr6_8_10Mbp_output.bam > Chr6_8_10Mbp_output.tbl
bgzip -c Chr6_8_10Mbp_output.tbl > Chr6_8_10Mbp_output.tbl.gz
python3.13 {SOFTWARE_LOCAL}/StainedGlass-0.6/workflow/scripts/refmt.py --window 2000 --fai Chr6_8_10Mbp.fa.fai --full Chr6_8_10Mbp_output.full.tbl.gz Chr6_8_10Mbp_output.tbl.gz Chr6_8_10Mbp_output.bed.gz
mkdir -p results/Chr6_8_10Mbp_figures/pdfs
mkdir -p results/Chr6_8_10Mbp_figures/pngs
Rscript {SOFTWARE_LOCAL}/StainedGlass-0.6/workflow/scripts/aln_plot.R -b Chr6_8_10Mbp_output.bed.gz --threads 88 --prefix Chr6_8_10Mbp


bedtools getfasta -fi ../Lotus_GifuT2T_v1.0.fasta -bed Chr6_25_28Mbp.bed > Chr6_25_28Mbp.fa
samtools faidx Chr6_25_28Mbp.fa
bedtools makewindows -g Chr6_25_28Mbp.fa.fai -w 2000 > Chr6_25_28Mbp_output.bed
bedtools getfasta -fi Chr6_25_28Mbp.fa -bed Chr6_25_28Mbp_output.bed > Chr6_25_28Mbp_output.2000.fasta

minimap2 -f 1000 -s 400 -ax ava-ont -d Chr6_25_28Mbp_output.fasta.mmi Chr6_25_28Mbp_output.2000.fasta
minimap2 -t 88 -f 10000 -s 400 -ax ava-ont --dual=yes --eqx Chr6_25_28Mbp_output.fasta.mmi Chr6_25_28Mbp_output.2000.fasta | samtools sort -m 4G -o Chr6_25_28Mbp_output.bam
python3.13 {SOFTWARE_LOCAL}/StainedGlass-0.6/workflow/scripts/samIdentity.py --threads 88 --matches 400 --header Chr6_25_28Mbp_output.bam > Chr6_25_28Mbp_output.tbl
bgzip -c Chr6_25_28Mbp_output.tbl > Chr6_25_28Mbp_output.tbl.gz
python3.13 {SOFTWARE_LOCAL}/StainedGlass-0.6/workflow/scripts/refmt.py --window 2000 --fai Chr6_25_28Mbp.fa.fai --full Chr6_25_28Mbp_output.full.tbl.gz Chr6_25_28Mbp_output.tbl.gz Chr6_25_28Mbp_output.bed.gz
mkdir -p results/Chr6_25_28Mbp_figures/pdfs
mkdir -p results/Chr6_25_28Mbp_figures/pngs
Rscript {SOFTWARE_LOCAL}/StainedGlass-0.6/workflow/scripts/aln_plot.R -b Chr6_25_28Mbp_output.bed.gz --threads 88 --prefix Chr6_25_28Mbp