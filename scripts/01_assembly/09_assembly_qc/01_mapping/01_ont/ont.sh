minimap2 -ax map-ont {PROJ_LOTUS_ZC}/Gifu/00.assemble_result/Gifu_v1.0.fasta {PROJ_LOTUS_ZC}/data/Gifu/ont/Gifu-jing-ye.pass.ul.fa -o ont.sam
samtools sort -@ 88 -o ont.bam -T ont_tmp.ali ont.sam
samtools flagstat ont.sam > ont.flagstat

samtools view -h -F 2304 ont.bam -o primary_only.bam

samtools view -b -q 30 -o ont_q30.bam ont.bam

bedtools coverage -a {PROJ_04_LOTUS_GENOME}/00_assembly_result/GifuT2T/Gifu_100k.bed -b ont_q30.bam -mean > gifu_ont_100k_coverage.txt
bedtools coverage -a {PROJ_04_LOTUS_GENOME}/00_assembly_result/GifuT2T/Gifu_10k.bed -b ont_q30.bam -mean > gifu_ont_10k_coverage.txt