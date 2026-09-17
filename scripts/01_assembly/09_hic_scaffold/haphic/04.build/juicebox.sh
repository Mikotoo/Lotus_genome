#!/bin/bash
#CSUB -J haphic
#CSUB -q c01
#CSUB -o haphic.out
#CSUB -e haphic.error
#CSUB -n 88
#CSUB -R span[hosts=1]

ln -s {PROJ_LOTUS_ZC}/Gifu/05.haphic/asm.fa .
samtools faidx asm.fa
{SOFTWARE_ZC}/HapHiC-main/scripts/../utils/juicer pre -a -q 1 -o out_JBAT {PROJ_LOTUS_ZC}/Gifu/05.haphic/HiC.filtered.bam scaffolds.raw.agp asm.fa.fai >out_JBAT.log 2>&1
(java -Djava.awt.headless=true -jar -Xmx32G {SOFTWARE_ZC}/HapHiC-main/scripts/../utils/juicer_tools.1.9.9_jcuda.0.8.jar pre out_JBAT.txt out_JBAT.hic.part <(cat out_JBAT.log | grep PRE_C_SIZE | awk '{print $2" "$3}')) && (mv out_JBAT.hic.part out_JBAT.hic)
