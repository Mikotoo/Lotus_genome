#!/bin/bash
#CSUB -J haphic
#CSUB -q c01
#CSUB -o haphic.out
#CSUB -e haphic.error
#CSUB -n 88
#CSUB -R span[hosts=1]

{SOFTWARE_ZC}/HapHiC-main/haphic pipeline asm.fa HiC.filtered.bam 6 --threads 48 --processes 48 --quick_view
#{SOFTWARE_ZC}/HapHiC-main/haphic pipeline asm.fa HiC.filtered.bam 6 --threads 48 --processes 48 --inflation 3