#!/bin/bash
# Reproduce panomap zymo multi-class classification from the Zenodo zymo/ tarball.
#
# Usage:
#   ./run_zymo.sh <extracted-zymo-dir> <panomap-binary> [out-dir]
#
# Env overrides: BUILD=pggb|mc  N_VALUES="1 2 4 8"  CHUNKS="1 2 4 8"  PY=python3
#
# Produces <out-dir>/results/summary/macro_summary.tsv (per-(N,chunks) macro_F1,
# mapped_frac) + per_species_metrics.tsv + confusion matrices.
set -euo pipefail

ZYMO=${1:?usage: run_zymo.sh <zymo-dir> <panomap-bin> [outdir]}
PANOMAP=${2:?usage: run_zymo.sh <zymo-dir> <panomap-bin> [outdir]}
OUT=${3:-./build/zymo}
PY=${PY:-python3}
EVAL="$(cd "$(dirname "$0")" && pwd)/eval/eval_zymo.py"

# ---- fixed zymo settings (paper Supplementary S4) ----
PORE=r9.4
IDX_T=16
MAP_T=64
MAP_ARGS="--seed-freq-cap p0.97 --map-sc-min-anchors 2"   # default path-chain

BUILD=${BUILD:-pggb}
N_VALUES=${N_VALUES:-1 2 4 8}
CHUNKS=${CHUNKS:-1 2 4 8}
ON="$ZYMO/reads/sigmoni-bacteria.blow5"
TRUTH="$ZYMO/reads/sigmoni-zymo.read-species.tsv"

mkdir -p "$OUT/idx"
for N in $N_VALUES; do
  GFA="$ZYMO/graphs/zymo-${BUILD}-${N}.gfa"
  IDX="$OUT/idx/zymo-${BUILD}-${N}.pirx"
  [ -f "$IDX" ] || "$PANOMAP" index "$GFA" -m "$PORE" -t "$IDX_T" -o "$IDX"
  for C in $CHUNKS; do
    # eval_zymo.py walks <out>/results/runs/zymo-N{N}/chunks_{C}/on_target.gaf
    CELL="$OUT/results/runs/zymo-N${N}/chunks_${C}"
    mkdir -p "$CELL"
    "$PANOMAP" map --index "$IDX" $MAP_ARGS \
      --max-chunks "$C" -t "$MAP_T" "$ON" -o "$CELL/on_target.gaf" 2> "$CELL/on.log"
  done
done

"$PY" "$EVAL" "$OUT" --truth "$TRUTH"
echo "=== macro_summary ($OUT/results/summary/macro_summary.tsv) ==="
column -t "$OUT/results/summary/macro_summary.tsv"
