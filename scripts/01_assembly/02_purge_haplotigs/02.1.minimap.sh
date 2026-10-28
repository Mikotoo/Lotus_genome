#2.1 把hifi reads回帖到1.1的初步组装结果上

minimap2 -ax map-pb Gifu.hifiasm.fasta {hifiReads}.fastq.gz \
|samtools view -hF 256 - \
|samtools sort -@ 20 -m 1G -o aligned.bam -T tmp.ali
