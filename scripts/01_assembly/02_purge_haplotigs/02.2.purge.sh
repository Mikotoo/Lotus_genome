#2.2 根据比对结果，去除单倍型冗余contig

purge_haplotigs readhist -b aligned.bam -g Gifu.hifiasm.fasta -t 88
purge_haplotigs contigcov -i aligned.bam.200.gencov -o coverage_stats.csv -l 15 -m 100 -h 185
purge_haplotigs purge -g Gifu.hifiasm.fasta -c coverage_stats.csv -b aligned.bam -t 88 -a 50

