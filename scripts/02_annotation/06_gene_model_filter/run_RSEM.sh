#!/bin/bash
#CSUB -J rsem
#CSUB -q c01
#CSUB -o rsem.out
#CSUB -e rsem.error
#CSUB -n 64
#CSUB -R span[host=1]


mkdir rsem_outdir
{CONDA_ZC}/envs/annotation/bin/align_and_estimate_abundance.pl --transcripts reference.fasta --seqType fq --left {PROJ_LOTUS_ZC}/Gifu/09.annotation/08.geta_filter/Gifu_RNA_all_R1.fastq.gz --right {PROJ_LOTUS_ZC}/Gifu/09.annotation/08.geta_filter/Gifu_RNA_all_R2.fastq.gz --est_method RSEM --aln_method bowtie --gene_trans_map gene.map --prep_reference --output_dir rsem_outdir --thread_count 64
samtools sort -@ 64 -o sorted.bam rsem_outdir/bowtie.bam
awk '{print $1"\t1\t"$2}' reference.fasta.fai > gene.bed
awk '{print $1"\t"$2}' reference.fasta.fai > transcriptome_chromSizes.txt
bedtools coverage -a gene.bed -b sorted.bam -sorted -g transcriptome_chromSizes.txt > cov.bed
