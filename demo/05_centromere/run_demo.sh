#!/bin/bash
# Demonstration of build_supercluster_unionfind.py on a three-supercluster example.
#
# input/    clusters.list and pairs.jaccard_ge0.7.tsv (hand-made, tiny)
# work/     scratch directory the script writes into
# expected_output/  the two files the script produces, kept for comparison
#
# Usage:  bash run_demo.sh
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
SCRIPT="$REPO/scripts/05_centromere/01_supercluster/blast/build_supercluster_unionfind.py"

rm -rf "$HERE/work"
mkdir -p "$HERE/work"
cp "$HERE/input/clusters.list" "$HERE/input/pairs.jaccard_ge0.7.tsv" "$HERE/work/"

python3 "$SCRIPT" --outdir "$HERE/work" --jac 0.7

echo
echo "--- cluster2super.tsv ---"
cat "$HERE/work/cluster2super.tsv"
echo
echo "--- super_summary.tsv ---"
cat "$HERE/work/super_summary.tsv"
echo
echo "Compare with $HERE/expected_output/"
