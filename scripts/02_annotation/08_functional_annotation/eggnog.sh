#!/bin/bash
#CSUB -J eggnog
#CSUB -q c01
#CSUB -o eggnog.out
#CSUB -e eggnog.error
#CSUB -n 64
#CSUB -R span[hosts=1]

python {SOFTWARE_ZC}/eggnog-mapper-2.1.13/emapper.py -m hmmer -d Eukaryota -i {PROJ_04_LOTUS_GENOME}/00_assembly_result/GifuT2T/Lotus_GifuT2T.pep.fa --output Lotus_Gifu --cpu 64