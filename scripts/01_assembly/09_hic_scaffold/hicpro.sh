#!/bin/bash
#CSUB -J hicpro
#CSUB -q c01
#CSUB -o hicpro.out
#CSUB -e hicpro.error
#CSUB -n 88
#CSUB -R span[hosts=1]

{SOFTWARE_ZC}/HiC-Pro/HiC-Pro_3.1.0/bin/HiC-Pro --input rawreads -c config-hicpro.txt --output process