#1.1 contig组装

hifiasm -t 64 -o Gifu.sam --telo-m CCCTAAA --ul {ontReads}.fq.gz {hifiReads}.fastq.gz

gfatools gfa2fa Gifu.sam.bp.p_ctg.gfa > Gifu.hifiasm.fasta

samtools faidx Gifu.hifiasm.fasta