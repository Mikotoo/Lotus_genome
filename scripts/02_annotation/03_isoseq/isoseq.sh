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

#> clustered_bams.txt

#for i in {1..4}; do
#    sample_var="sample${i}"
#    primer_var="primer${i}"
#    sample=$(eval echo \$$sample_var)
#    primer=$(eval echo \$$primer_var)

#    prefix=$(basename "$sample" .hifi_reads.bam)

#    echo "Processing $prefix"
#    isoseq3 refine "$sample" "$primer" "${prefix}.refined.bam"
#    isoseq3 cluster "${prefix}.refined.bam" "${prefix}.clustered.bam" --verbose --use-qvs
#    echo "${prefix}.clustered.bam" >> clustered_bams.txt
#done

#bamtools merge $(cat clustered_bams.txt | awk '{print "-in " $1}') -out all_samples.clustered.bam

samtools fastq Gifu.clustered.bam  > Gifu.clustered.fastq
minimap2 -ax splice -uf -k14 {PROJ_LOTUS_ZC}/Gifu/00.assemble_result/Gifu_v1.0.fasta Gifu.clustered.fastq > Gifu_isoforms.sam
samtools sort -o Gifu_isoforms.sorted.bam Gifu_isoforms.sam
samtools index Gifu_isoforms.sorted.bam

#python {SOFTWARE_ZC}/tama/tama_collapse_py3.py -s isoforms.sorted.bam -b BAM -f updated_reference.fa -p isoseq
#python {SOFTWARE_ZC}/tama/tama_go/format_converter/tama_convert_bed_gtf_ensembl_no_cds.py isoseq.bed isoseq.gtf

#gffread isoseq.gtf -g {HOME_ZC}/Project/GZ/annotation/liftoff/updated_reference.fa -w transcripts.fa

## << TransDecoder >>

#{SOFTWARE_ZC}/TransDecoder/util/gtf_genome_to_cdna_fasta.pl isoseq.gtf {HOME_ZC}/Project/GZ/annotation/liftoff/updated_reference.fa > isoseq.cdna.fasta
#{SOFTWARE_ZC}/TransDecoder/util/gtf_to_alignment_gff3.pl isoseq.gtf > isoseq.gff3


#TransDecoder.LongOrfs -t transcripts.fa
#TransDecoder.Predict -t transcripts.fa


#{SOFTWARE_ZC}/TransDecoder/util/cdna_alignment_orf_to_genome_orf.pl \
#transcripts.fa.transdecoder.gff3 isoseq.gff3 isoseq.cdna.fasta > isoseq.genome_orf.gff3