#!/bin/bash
# Reproduce panomap yeast binary classification from the Zenodo yeast/ tarball.
#
# Usage:
#   ./run_yeast.sh <extracted-yeast-dir> <panomap-binary> [out-dir]
#
# Env overrides: BUILD=pggb|mc  N_VALUES="1 2 4 8"  CHUNKS="1 2 4 8"  PY=python3
#
# Produces <out-dir>/results.tsv with per-(N,chunks) TP/FP/FN/TN + P/R/spec/F1.
set -euo pipefail

YEAST=${1:?usage: run_yeast.sh <yeast-dir> <panomap-bin> [outdir]}
PANOMAP=${2:?usage: run_yeast.sh <yeast-dir> <panomap-bin> [outdir]}
OUT=${3:-./build/yeast}
PY=${PY:-python3}
EVAL="$(cd "$(dirname "$0")" && pwd)/eval/eval_binary.py"

# ---- fixed yeast settings (paper Supplementary S4) ----
PORE=r9.4
IDX_T=16
MAP_T=64
CHAINER=path-chain
PATHCHAIN="--chain-pen-gap 1.4 --chain-pen-skip 0.1"
GATE=""                       # yeast uses single-chain gate defaults
ON_TOTAL=40000
OFF_TOTAL=40000

BUILD=${BUILD:-pggb}
N_VALUES=${N_VALUES:-1 2 4 8}
CHUNKS=${CHUNKS:-1 2 4 8}
ON="$YEAST/reads/yeast-r9-40k.blow5"
OFF="$YEAST/reads/ecoli-sim-r9-4khz-40k-8kb.blow5"

mkdir -p "$OUT/idx"
"$PY" "$EVAL" --header > "$OUT/results.tsv"

for N in $N_VALUES; do
  GFA="$YEAST/graphs/yeast-${BUILD}-${N}.gfa"
  IDX="$OUT/idx/yeast-${BUILD}-${N}.pirx"
  [ -f "$IDX" ] || "$PANOMAP" index "$GFA" -m "$PORE" -t "$IDX_T" -o "$IDX"
  for C in $CHUNKS; do
    CELL="$OUT/runs/yeast-${BUILD}-N${N}/chunks_${C}"
    mkdir -p "$CELL"
    "$PANOMAP" map --index "$IDX" --chainer "$CHAINER" $PATHCHAIN $GATE \
      --max-chunks "$C" -t "$MAP_T" "$ON"  -o "$CELL/on_target.gaf"  2> "$CELL/on.log"
    "$PANOMAP" map --index "$IDX" --chainer "$CHAINER" $PATHCHAIN $GATE \
      --max-chunks "$C" -t "$MAP_T" "$OFF" -o "$CELL/off_target.gaf" 2> "$CELL/off.log"
    "$PY" "$EVAL" --on "$CELL/on_target.gaf" --off "$CELL/off_target.gaf" \
      --on-total "$ON_TOTAL" --off-total "$OFF_TOTAL" \
      --name "yeast-${BUILD}-N${N}" --chunks "$C" >> "$OUT/results.tsv"
  done
done

echo "=== results ($OUT/results.tsv) ==="
column -t "$OUT/results.tsv"
