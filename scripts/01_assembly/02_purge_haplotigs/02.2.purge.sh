#!/bin/bash
#CSUB -J purge
#CSUB -q c01
#CSUB -o purge.out
#CSUB -e purge.error
#CSUB -n 88
#CSUB -R span[hosts=1]

#purge_haplotigs readhist -b aligned.bam -g {PROJ_LOTUS_ZC}/Gifu/01.hifiasm/Gifu.hifiasm.fasta -t 88
#purge_haplotigs contigcov -i aligned.bam.200.gencov -o coverage_stats.csv -l 15 -m 100 -h 185
purge_haplotigs purge -g {PROJ_LOTUS_ZC}/Gifu/01.hifiasm/Gifu.hifiasm.fasta -c coverage_stats.csv -b aligned.bam -t 88 -a 50

