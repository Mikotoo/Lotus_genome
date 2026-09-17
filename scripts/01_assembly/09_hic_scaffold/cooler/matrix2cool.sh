#!/bin/bash
#CSUB -J cooler
#CSUB -q c01
#CSUB -o cooler.out
#CSUB -e cooler.error
#CSUB -n 64
#CSUB -R span[hosts=1]

cat {PROJ_04_LOTUS_GENOME}/02_assembly/Gifu/09.hic/process/hic_results/matrix/Gifu/raw/200000/Gifu_200000.matrix \
| awk '$1<=$2{print $0}' > {PROJ_04_LOTUS_GENOME}/02_assembly/Gifu/09.hic/process/hic_results/matrix/Gifu/raw/200000/Gifu_200000_up.matrix

cooler load -f coo --one-based --count-as-float \
{PROJ_04_LOTUS_GENOME}/02_assembly/Gifu/09.hic/process/hic_results/matrix/Gifu/raw/200000/Gifu_200000_abs.bed \
{PROJ_04_LOTUS_GENOME}/02_assembly/Gifu/09.hic/process/hic_results/matrix/Gifu/raw/200000/Gifu_200000_up.matrix \
Gifu_200kb.cool


cat {PROJ_04_LOTUS_GENOME}/02_assembly/Gifu/09.hic/process/hic_results/matrix/Gifu/raw/100000/Gifu_100000.matrix \
| awk '$1<=$2{print $0}' > {PROJ_04_LOTUS_GENOME}/02_assembly/Gifu/09.hic/process/hic_results/matrix/Gifu/raw/100000/Gifu_100000_up.matrix

cooler load -f coo --one-based --count-as-float \
{PROJ_04_LOTUS_GENOME}/02_assembly/Gifu/09.hic/process/hic_results/matrix/Gifu/raw/100000/Gifu_100000_abs.bed \
{PROJ_04_LOTUS_GENOME}/02_assembly/Gifu/09.hic/process/hic_results/matrix/Gifu/raw/100000/Gifu_100000_up.matrix \
Gifu_100kb.cool


cat {PROJ_04_LOTUS_GENOME}/02_assembly/Gifu/09.hic/process/hic_results/matrix/Gifu/raw/40000/Gifu_40000.matrix \
| awk '$1<=$2{print $0}' > {PROJ_04_LOTUS_GENOME}/02_assembly/Gifu/09.hic/process/hic_results/matrix/Gifu/raw/40000/Gifu_40000_up.matrix

cooler load -f coo --one-based --count-as-float \
{PROJ_04_LOTUS_GENOME}/02_assembly/Gifu/09.hic/process/hic_results/matrix/Gifu/raw/40000/Gifu_40000_abs.bed \
{PROJ_04_LOTUS_GENOME}/02_assembly/Gifu/09.hic/process/hic_results/matrix/Gifu/raw/40000/Gifu_40000_up.matrix \
Gifu_40kb.cool



cooler create --mode mcool -o Gifu.mcool
cooler cp Gifu_200kb.cool Gifu.mcool::resolutions/200000
cooler cp Gifu_100kb.cool Gifu.mcool::resolutions/100000
cooler cp Gifu_40kb.cool Gifu.mcool::resolutions/40000


for resolution_uri in $(cooler ls Gifu.mcool); do
    cooler balance \
        --ignore-diags 2 \
        --force \
        "$resolution_uri"
done