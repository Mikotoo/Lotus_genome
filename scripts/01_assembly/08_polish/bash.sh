#8.1 分别用ont，hifi reads polish 2轮

polish.sh 88 2 Gifu_v0.9.fa {ontReads}.fq.gz ccs.merylDB ontPolish
polish.sh 88 2 ontPolish.iter_2.consensus.fasta {hifiReads}.fq.gz ccs.merylDB hifiPolish
