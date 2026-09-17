#!/bin/bash
#CSUB -J telomere
#CSUB -q c01
#CSUB -o run_telomere.out
#CSUB -e run_telomere.error
#CSUB -n 32
#CSUB -R span[hosts=1]
#
# Gifu telomere analysis — single entry point.
#
#   csub run_telomere.sh                    # full pipeline (compute + one plot pass)
#   TELOMERE_STEPS=plot bash run_telomere.sh   # re-plot only
#   TELOMERE_STEPS=compute bash run_telomere.sh  # find-tel + coverage + alignments
#
# Python: telomere.py  (subcommands: find-tel | coverage | alignments | plot | all)

set -euo pipefail

PROJECT_ROOT="{PROJ_04_LOTUS_GENOME}"
# csub copies the script to ~/.cbsched; BASH_SOURCE then points there, not at telomere.py.
SCRIPT_DIR="${PROJECT_ROOT}/03_annotation/Gifu/01.tel"
WORKDIR="${TELOMERE_WORKDIR:-${SCRIPT_DIR}}"
PREFIX="${TELOMERE_PREFIX:-GifuT2T}"
STEPS="${TELOMERE_STEPS:-all}"

PYTHON="${TELOMERE_PYTHON:-{CONDA_PY3}/bin/python3}"
SAMTOOLS="${TELOMERE_SAMTOOLS:-{SOFTWARE_ZC}/TeloComp-1.0.0/Dependencies/samtools-1.18/samtools}"

FASTA_ASM="${TELOMERE_FASTA_ASM:-${PROJECT_ROOT}/00_assembly_result/GifuT2T/Lotus_GifuT2T_v1.0.fasta}"
FASTA_PLOT="${TELOMERE_FASTA_PLOT:-${PROJECT_ROOT}/00_assembly_result/GifuT2T/Lotus_GifuT2T_v1.0.fasta}"
ONT_BAM="${TELOMERE_ONT_BAM:-${PROJECT_ROOT}/02_assembly/Gifu/08.quality/01.mapping/01.ont/ont_q30.bam}"
HIFI_BAM="${TELOMERE_HIFI_BAM:-${PROJECT_ROOT}/02_assembly/Gifu/08.quality/01.mapping/02.hifi/ccs_q30.bam}"
BEDTOOLS="${TELOMERE_BEDTOOLS:-{CONDA_PY3}/bin/bedtools}"

BIN_SIZE="${TELOMERE_BIN_SIZE:-100}"
PLOT_FLANK="${TELOMERE_PLOT_FLANK:-10000}"
THREADS="${TELOMERE_THREADS:-16}"

BASE_ARGS=(
  --workdir "${WORKDIR}"
  --prefix "${PREFIX}"
  --bin-size "${BIN_SIZE}"
  --plot-flank "${PLOT_FLANK}"
  --fasta-asm "${FASTA_ASM}"
  --fasta-plot "${FASTA_PLOT}"
  --ont-bam "${ONT_BAM}"
  --hifi-bam "${HIFI_BAM}"
  --bedtools "${BEDTOOLS}"
)

index_bam() {
  local bam="$1"
  if [[ ! -f "${bam}.bai" ]]; then
    echo "Indexing ${bam} ..."
    "${SAMTOOLS}" index -@"${THREADS}" "${bam}"
  fi
}

cd "${WORKDIR}"
echo "WORKDIR : ${WORKDIR}"
echo "PREFIX  : ${PREFIX}"
echo "STEPS   : ${STEPS}"

case "${STEPS,,}" in
  all)
    if [[ ! -f "${FASTA_ASM}.fai" ]]; then
      echo "Building FAI: ${FASTA_ASM}.fai"
      "${SAMTOOLS}" faidx "${FASTA_ASM}"
    fi
    index_bam "${ONT_BAM}"
    index_bam "${HIFI_BAM}"
    "${PYTHON}" "${SCRIPT_DIR}/telomere.py" all "${BASE_ARGS[@]}" --reuse-units
    ;;
  compute)
    index_bam "${ONT_BAM}"
    index_bam "${HIFI_BAM}"
    for step in find-tel coverage alignments; do
      "${PYTHON}" "${SCRIPT_DIR}/telomere.py" "${step}" "${BASE_ARGS[@]}" --reuse-units
    done
    ;;
  plot)
    "${PYTHON}" "${SCRIPT_DIR}/telomere.py" plot "${BASE_ARGS[@]}"
    ;;
  find-tel|find_tel|coverage|alignments|align|reads)
    step="${STEPS,,}"
    step="${step//_/-}"
    [[ "${step}" == "find_tel" ]] && step="find-tel"
    [[ "${step}" == "align" || "${step}" == "reads" ]] && step="alignments"
    if [[ "${step}" != "find-tel" ]]; then
      index_bam "${ONT_BAM}"
      index_bam "${HIFI_BAM}"
    fi
    extra=()
    [[ "${step}" == "find-tel" ]] && extra=(--reuse-units)
    "${PYTHON}" "${SCRIPT_DIR}/telomere.py" "${step}" "${BASE_ARGS[@]}" "${extra[@]}"
    ;;
  *)
    echo "ERROR: unknown TELOMERE_STEPS=${STEPS}" >&2
    echo "  Use: all | compute | plot | find-tel | coverage | alignments" >&2
    exit 1
    ;;
esac

echo "Done. Tables: ${WORKDIR}/${PREFIX}.*  Figures: ${WORKDIR}/plots/"
