#!/bin/bash
#CSUB -J nucmer
#CSUB -q c01
#CSUB -o nucmer.out
#CSUB -e nucmer.error
#CSUB -n 64
#CSUB -R span[hosts=1]

#nucmer --prefix ref {PROJ_LOTUS_ZC}/Gifu/05_3.mummer/Gifu_ref.fa {PROJ_LOTUS_ZC}/MG20/05_3.mummer/MG20_ref.fa
#delta-filter -i 89 -l 1000 -1 ref.delta > ref.filter.delta
#show-coords -THrd ref.filter.delta > ref.filter.delta.coords

mummerplot ref.filter.delta -R {PROJ_LOTUS_ZC}/Gifu/05_3.mummer/Gifu_ref.fa -Q {PROJ_LOTUS_ZC}/MG20/05_3.mummer/MG20_ref.fa --layout --png --large