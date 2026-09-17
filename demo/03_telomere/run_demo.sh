#!/bin/bash
# Demonstration of the telomere pipeline on a short synthetic chromosome.
#
# input/demo_chr.fa        2 kb sequence: 40 x TTTAGGG, filler, 40 x CCCTAAA
# input/demo_chr.fa.fai    FASTA index required by telomere.py
# expected_output/         the four files the run produces
#
# Only the find-tel subcommand is exercised; coverage, alignments and plot need BAM
# files and a full assembly, which are not part of a small demo.
#
# Usage:  bash run_demo.sh
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
SCRIPT="$REPO/scripts/03_telomere/telomere.py"

rm -rf "$HERE/work"
mkdir -p "$HERE/work"

python3 "$SCRIPT" find-tel \
  --fasta-asm "$HERE/input/demo_chr.fa" \
  --workdir "$HERE/work" \
  --prefix demo

echo
echo "--- demo.telomere_coords.tsv ---"
cat "$HERE/work/demo.telomere_coords.tsv"
echo
echo "Compare with $HERE/expected_output/"
