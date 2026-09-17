#!/bin/bash
#CSUB -J count
#CSUB -q c01
#CSUB -o count.out
#CSUB -e count.error
#CSUB -n 64
#CSUB -R span[hosts=1]