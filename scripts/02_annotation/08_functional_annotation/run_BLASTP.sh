#!/bin/bash
cd $PBS_O_WORKDIR
ParaFly -c cmd.list -CPU 0
cat blast*.out.tmp > BLASTP.OUT.TMP
rm blast*.out.tmp
