#!/bin/bash
#CSUB -J blast
#CSUB -q c01
#CSUB -o blast.out
#CSUB -e blast.error
#CSUB -n 88
#CSUB -R span[hosts=1]

makeblastdb -in douke_pep.fasta -dbtype prot -out douke_pep.fasta