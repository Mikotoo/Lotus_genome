#!/bin/bash
#CSUB -J blast
#CSUB -q c02
#CSUB -o blast.out
#CSUB -e blast.error
#CSUB -n 64
#CSUB -R span[hosts=1]



GENOME={PROJ_04_LOTUS_GENOME}/00_assembly_result/GifuT2T/Lotus_GifuT2T_v1.0.fasta
DB={PROJ_04_LOTUS_GENOME}/03_annotation/Gifu/02.cent/blast/Gifu_genome_db

#makeblastdb -in "$GENOME" -dbtype nucl -out "$DB"


#blastn \
#  -task blastn \
#  -db "$DB" \
#  -query Gifu_cluster_reps.fa \
#  -evalue 1e-10 \
#  -perc_identity 80 \
#  -dust no \
# -outfmt "6 qseqid sseqid sstart send pident length bitscore" \
#  -num_threads 64 \
# > Gifu_cluster_reps.blast


bash cluster_density.sh Gifu_cluster_reps.blast {PROJ_04_LOTUS_GENOME}/00_assembly_result/GifuT2T/Lotus_GifuT2T_v1.0.fasta supercluster
