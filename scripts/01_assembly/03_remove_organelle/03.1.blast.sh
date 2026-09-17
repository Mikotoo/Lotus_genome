#!/bin/bash
#CSUB -J blastn
#CSUB -q c01
#CSUB -o blastn.out
#CSUB -e blastn.error
#CSUB -n 88
#CSUB -R span[hosts=1]


#1 提取最短的几条序列，比对到nt库，输出fmt5格式
blastn -query contigs.fa -out nt.blast.xml -db {SOFTWARE_ZC}/nt_lib/nt -outfmt 5 -evalue 1e-5 -num_threads 88

python 03.2.get_gi.py
python 03.3.get_name.py