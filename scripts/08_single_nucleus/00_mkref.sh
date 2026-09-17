#!/bin/bash
#CSUB -J mkref
#CSUB -q c01
#CSUB -o mkref.out
#CSUB -e mkref.error
#CSUB -n 64
#CSUB -R span[hosts=1]

reffasta={PROJ_LOTUS_HX}/00_assembly_result/GifuT2T/Lotus_GifuT2T_v1.0.fasta
refgtf={PROJ_LOTUS_HX}/00_assembly_result/GifuT2T/Lotus_GifuT2T.gtf

cellranger mkref \
  --genome=Gifu \
  --fasta=$reffasta \
  --genes=$refgtf \
  --output-dir=./00_ref \
  --nthreads=64