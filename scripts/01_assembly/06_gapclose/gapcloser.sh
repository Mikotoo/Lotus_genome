#!/bin/bash
#CSUB -J tgsgapcloser
#CSUB -q c01
#CSUB -o tgsgapcloser.out
#CSUB -e tgsgapcloser.error
#CSUB -n 64
#CSUB -R span[host=1]

{SOFTWARE_ZC}/TGS-GapCloser2/tgsgapcloser2 --scaff chr2.fa \
--reads {PROJ_LOTUS_ZC}/data/Gifu/ont/Gifu-jing-ye.pass.ul.fa \
--output ont_fill_out \
--ne --thread 64
