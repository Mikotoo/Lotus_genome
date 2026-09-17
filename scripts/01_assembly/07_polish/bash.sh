#!/bin/bash
#CSUB -J polish
#CSUB -q c01
#CSUB -o polish.out
#CSUB -e polish.error
#CSUB -n 88
#CSUB -R span[host=1]

#{PROJ_LOTUS_ZC}/Gifu/07.polish/polish.sh 88 2 Gifu_v0.9.fa {PROJ_LOTUS_ZC}/data/Gifu/ont/Gifu-jing-ye.pass.ul.fa {PROJ_LOTUS_ZC}/Gifu/07.polish/ccs.merylDB ontPolish
{PROJ_LOTUS_ZC}/Gifu/07.polish/polish.sh 10 2 ontPolish.iter_2.consensus.fasta {PROJ_LOTUS_ZC}/data/Gifu/hifi/Gifu_hifi.fastq.gz {PROJ_LOTUS_ZC}/Gifu/07.polish/ccs.merylDB hifiPolish
