#!/bin/bash
#CSUB -J EDTA
#CSUB -q c01
#CSUB -o EDTA.out
#CSUB -e EDTA.error
#CSUB -n 264
#CSUB -R span[ptile=88]

EDTA.pl --genome Gifu_v1.0.fasta \
--overwrite 1 --sensitive 1 --anno 1 --threads 264 --force 1