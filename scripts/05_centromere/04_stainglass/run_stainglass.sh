#!/bin/bash
#CSUB -J process
#CSUB -q c01
#CSUB -o process.out
#CSUB -e process.error
#CSUB -n 88
#CSUB -R span[hosts=1]

# StainedGlass 0.6 self-identity run for the six Gifu centromere regions.
#
# This is the loop form of the six per-region scripts that were run one
# directory at a time (Gifu_Chr1 … Gifu_Chr6). Each directory holds the region
# FASTA produced by get.sh: <region>.fa, <region>.fa.fai and <region>_2k.fasta.
#
# Run from 04_stainglass/, after get.sh:
#   bash run_stainglass.sh

SG=~/software/StainedGlass-0.6/workflow/scripts

for dir in Gifu_Chr*; do
    [ -d "$dir" ] || continue
    region="$dir"
    echo ">>> ${region}"
    cd "$dir"

    minimap2 -f 1000 -s 400 -ax ava-ont -d ${region}_2k.mmi ${region}_2k.fasta
    minimap2 -t 88 -f 10000 -s 400 -ax ava-ont --dual=yes --eqx ${region}_2k.mmi ${region}_2k.fasta | samtools sort -m 4G -o ${region}.bam
    python ${SG}/samIdentity.py --threads 88 --matches 400 --header ${region}.bam > ${region}.tbl
    bgzip -c ${region}.tbl > ${region}.tbl.gz
    python ${SG}/refmt.py --window 2000 --fai ${region}.fa.fai --full ${region}.full.tbl.gz ${region}.tbl.gz ${region}.bed.gz
    mkdir -p results/${region}_figures/pdfs
    mkdir -p results/${region}_figures/pngs
    Rscript ${SG}/aln_plot.R -b ${region}.bed.gz --threads 88 --prefix ${region}

    cd ..
done
