#!/bin/bash
# Reproduce panomap covid binary classification from the Zenodo covid/ tarball.
#
# Usage:
#   ./run_covid.sh <extracted-covid-dir> <panomap-binary> [out-dir]
#
# Env overrides (optional):
#   BUILD=pggb|mc     graph builder (default pggb)
#   N_VALUES="1 2 4 8"  pangenome sizes to run (default all)
#   CHUNKS="1 2 4 8"    max-chunks grid (default all)
#   PY=python3          python for the eval script
#
# Produces <out-dir>/results.tsv with per-(N,chunks) TP/FP/FN/TN + P/R/spec/F1.
set -euo pipefail

COVID=${1:?usage: run_covid.sh <covid-dir> <panomap-bin> [outdir]}
PANOMAP=${2:?usage: run_covid.sh <covid-dir> <panomap-bin> [outdir]}
OUT=${3:-./build/covid}
PY=${PY:-python3}
EVAL="$(cd "$(dirname "$0")" && pwd)/eval/eval_binary.py"

# ---- fixed covid settings (paper Supplementary S4) ----
PORE=r9.4
IDX_T=16
MAP_T=64
CHAINER=path-chain
PATHCHAIN="--chain-bw 100 --chain-max-dist 500 --chain-min-score 10 --chain-pen-gap 1.2 --chain-pen-skip 0.3"
GATE="--map-sc-min-anchors 2 --map-sc-ratio-lo 0.5 --map-sc-ratio-hi 1.7"
ON_TOTAL=40000
OFF_TOTAL=40000

BUILD=${BUILD:-pggb}
N_VALUES=${N_VALUES:-1 2 4 8}
CHUNKS=${CHUNKS:-1 2 4 8}
ON="$COVID/reads/covid-interartic-r9-40k.blow5"
OFF="$COVID/reads/ecoli-sim-r9-4khz-40k-8kb.blow5"

mkdir -p "$OUT/idx"
"$PY" "$EVAL" --header > "$OUT/results.tsv"

for N in $N_VALUES; do
  GFA="$COVID/graphs/covid-${BUILD}-${N}.gfa"
  IDX="$OUT/idx/covid-${BUILD}-${N}.pirx"
  [ -f "$IDX" ] || "$PANOMAP" index "$GFA" -m "$PORE" -t "$IDX_T" -o "$IDX"
  for C in $CHUNKS; do
    CELL="$OUT/runs/covid-${BUILD}-N${N}/chunks_${C}"
    mkdir -p "$CELL"
    "$PANOMAP" map --index "$IDX" --chainer "$CHAINER" $PATHCHAIN $GATE \
      --max-chunks "$C" -t "$MAP_T" "$ON"  -o "$CELL/on_target.gaf"  2> "$CELL/on.log"
    "$PANOMAP" map --index "$IDX" --chainer "$CHAINER" $PATHCHAIN $GATE \
      --max-chunks "$C" -t "$MAP_T" "$OFF" -o "$CELL/off_target.gaf" 2> "$CELL/off.log"
    "$PY" "$EVAL" --on "$CELL/on_target.gaf" --off "$CELL/off_target.gaf" \
      --on-total "$ON_TOTAL" --off-total "$OFF_TOTAL" \
      --name "covid-${BUILD}-N${N}" --chunks "$C" >> "$OUT/results.tsv"
  done
done

echo "=== results ($OUT/results.tsv) ==="
column -t "$OUT/results.tsv"
