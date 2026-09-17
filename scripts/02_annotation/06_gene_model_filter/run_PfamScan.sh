#!/bin/bash
#CSUB -J run_pfam
#CSUB -q c01
#CSUB -o run_pfam.out
#CSUB -e run_pfam.error
#CSUB -n 88
#CSUB -R span[host=1]

perl {SOFTWARE_ZC}/PfamScan/pfam_scan.pl -fasta  {PROJ_LOTUS_ZC}/Gifu/09.annotation/07.geta/Gifu_T2T.protein.fasta -dir {SOFTWARE_ZC}/PfamScan/data -cpu 88 -outfile pfam.out
