#!/bin/bash
#CSUB -J quarTeT
#CSUB -q c01
#CSUB -o quarTeT.out
#CSUB -e quarTeT.error
#CSUB -n 88
#CSUB -R span[host=1]

python ~/software/quarTeT/quartet.py GapFiller \
-d chr2.fa \
-g {PROJ_LOTUS_ZC}/Gifu_hifionly/04.tel_check/contigs.fa \
-p Gifu_gapfill \
-f 20000 -i 90 \
-t 88