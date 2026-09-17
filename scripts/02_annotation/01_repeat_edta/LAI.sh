#!/bin/bash
#CSUB -J LAI
#CSUB -q c01
#CSUB -o LAI.out
#CSUB -e LAI.error
#CSUB -n 88
#CSUB -R span[host=1]


LTR_retriever \
  -genome ChiHei_v1.0.fasta \
  -inharvest "ChiHei_v1.0.fasta.mod.harvest.combine.scn" \
  -infinder  "ChiHei_v1.0.fasta.mod.finder.combine.scn" \
  -threads 88