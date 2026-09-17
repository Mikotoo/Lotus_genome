#!/bin/bash
#CSUB -J blastn
#CSUB -q c01
#CSUB -o blastn.out
#CSUB -e blastn.error
#CSUB -n 88
#CSUB -R span[hosts=1]


blastn -query MG20_subsample_100k.fa -out MG20_subsample.blast.xml -db {SOFTWARE_ZC}/nt_lib/nt -outfmt 5 -evalue 1e-5 -num_threads 88

python blast_xml_species_summary.py \
-i MG20_subsample.blast.xml \
-a {HOME_ZC}/Project/Stero-seq/unmapped_seq_check/cutted_nucl_gb.accession2taxid \
-n {HOME_ZC}/Project/Stero-seq/unmapped_seq_check/names.dmp \
-o MG20_subsample

python 03.2.get_gi.py