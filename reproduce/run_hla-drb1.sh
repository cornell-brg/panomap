#!/bin/bash
# Reproduce panomap HLA-DRB1 leave-one-out recall from the Zenodo hla-drb1/ tarball.
#
# Usage:
#   ./run_hla-drb1.sh <extracted-hla-drb1-dir> <panomap-binary> <tier> [out-dir]
#     tier = 15 (near) | 12 (moderate) | 07 (far)
#
# Env overrides: CHUNKS="1 2 4 8"  PY=python3
#
# For the tier, indexes N=1 (GRCh38 linear) + N=48 + N=max train graphs,
# maps the tier's simulated reads against each, and reports per-(N,chunks) recall.
set -euo pipefail

DRB1=${1:?usage: run_hla-drb1.sh <hla-drb1-dir> <panomap-bin> <tier> [outdir]}
PANOMAP=${2:?usage: run_hla-drb1.sh <hla-drb1-dir> <panomap-bin> <tier> [outdir]}
G=${3:?tier: 15 | 12 | 07}
OUT=${4:-./build/hla-drb1/G$G}
PY=${PY:-python3}
EVAL="$(cd "$(dirname "$0")" && pwd)/eval/eval_drb1_recall.py"

# ---- fixed drb1 settings (paper Supplementary S4: defaults) ----
PORE=r9.4
IDX_T=16
MAP_T=64
CHAINER=path-chain
CHUNKS=${CHUNKS:-1 2 4 8}
READS="$DRB1/reads/drb1-G${G}.blow5"

# N=1 = shared GRCh38 linear; N=48 and N=max = per-tier train graphs.
# N values are derived from the graph filenames (drb1-G<tier>-train-<N>.gfa).
declare -A GFA
GFA[1]="$DRB1/graphs/drb1-grch38-linear.gfa"
for g in "$DRB1"/graphs/drb1-G${G}-train-*.gfa; do
  n=$(basename "$g" | sed -E 's/.*train-([0-9]+)\.gfa/\1/')
  GFA[$n]="$g"
done

mkdir -p "$OUT/idx"
for N in $(printf '%s\n' "${!GFA[@]}" | sort -n); do
  IDX="$OUT/idx/drb1-N${N}.pirx"
  [ -f "$IDX" ] || "$PANOMAP" index "${GFA[$N]}" -m "$PORE" -t "$IDX_T" -o "$IDX"
  for C in $CHUNKS; do
    CELL="$OUT/runs/N${N}/chunks_${C}"
    mkdir -p "$CELL"
    "$PANOMAP" map --index "$IDX" --chainer "$CHAINER" \
      --max-chunks "$C" -t "$MAP_T" "$READS" -o "$CELL/on_target.gaf" 2> "$CELL/on.log"
  done
done

echo "=== recall (tier G$G) ==="
"$PY" "$EVAL" "$OUT/runs"
