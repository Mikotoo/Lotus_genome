#!/bin/bash
#CSUB -J nucmer
#CSUB -q c01
#CSUB -o nucmer.out
#CSUB -e nucmer.error
#CSUB -n 64
#CSUB -R span[hosts=1]

nucmer --prefix Gifu_scaffold Gifu_ref.fa Gifu_contig.fa
delta-filter -i 89 -l 1000 -1 Gifu_scaffold.delta > Gifu_scaffold.filter.delta
show-coords -THrd Gifu_scaffold.filter.delta > Gifu_scaffold.filter.delta.coords

show-coords -THrd Gifu_scaffold.delta > Gifu_scaffold.delta.coords

mummerplot Gifu_scaffold.filter.delta -R Gifu_ref.fa -Q Gifu_contig.fa --layout --png --large