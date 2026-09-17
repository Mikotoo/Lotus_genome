#!/bin/bash
#CSUB -J cat
#CSUB -q c01
#CSUB -o cat.out
#CSUB -e cat.error
#CSUB -n 88
#CSUB -R span[host=1]

ParaFly -c cmd.list -CPU 80
cat blast*.out.tmp > BLASTP.OUT.TMP
rm blast*.out.tmp
