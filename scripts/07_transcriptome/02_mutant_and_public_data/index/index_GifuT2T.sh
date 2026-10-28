#!/bin/bash
#CSUB -J hisat_build
#CSUB -q c01
#CSUB -o hisat_build.out
#CSUB -e hisat_build.error
#CSUB -n 32
#CSUB -R span[hosts=1]

# 构建 Lotus japonicus GifuT2T 的 hisat2 索引

set -euo pipefail
export PATH={CONDA_PY3}/bin:$PATH

FASTA="{PROJ_04_LOTUS_GENOME}/00_assembly_result/GifuT2T/Lotus_GifuT2T_v1.0.fasta"
INDEX_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="${INDEX_DIR}/GifuT2T"

hisat2-build -p 32 "${FASTA}" "${OUT}"
echo "[DONE] index built at ${OUT}"
