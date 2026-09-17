#!/bin/bash
cd $PBS_O_WORKDIR
perl {SOFTWARE_ZC}/PfamScan/pfam_scan.pl -fasta  {PROJ_LOTUS_ZC}/Gifu/09.annotation/07.geta/Gifu_T2T.protein.fasta -dir {SOFTWARE_ZC}/Pfam -cpu 1 -outfile pfam.out
