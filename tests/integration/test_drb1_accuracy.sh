#!/bin/bash
# Integration test: DRB1 pangenome mapping accuracy
# Simulates 1000 reads, indexes DRB1 graph, maps, checks accuracy >= 90%
#
# Requires: panomap binary, squigulator, eval_canonical.py, benchmark data

set -euo pipefail

PANOMAP="${1:?Usage: $0 <panomap_binary> <workspace_root>}"
WORKSPACE="${2:?Usage: $0 <panomap_binary> <workspace_root>}"

SCRIPTS="$WORKSPACE/repo/scripts"
DATA="$WORKSPACE/data/benchmark/hla-drb1"
GFA="$DATA/drb1-pggb.gfa"
FASTA="$DATA/DRB1-3123.fa"
SQUIGULATOR="$WORKSPACE/infra/tools/squigulator"
PYTHON="${WORKSPACE}/infra/env/panomap-py/bin/python3"

NUM_READS=1000
READ_LEN=12000
SIM_SEED=42
MIN_ACCURACY=90

# Check inputs
for f in "$GFA" "$FASTA" "$SCRIPTS/eval_canonical.py"; do
  if [ ! -f "$f" ]; then
    echo "SKIP: missing $f"
    exit 0
  fi
done

if [ ! -x "$SQUIGULATOR" ]; then
  echo "SKIP: squigulator not found at $SQUIGULATOR"
  exit 0
fi

if [ ! -x "$PYTHON" ]; then
  echo "SKIP: python not found at $PYTHON"
  exit 0
fi

TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# 1. Simulate reads
"$SQUIGULATOR" -x dna-r10-prom --seed $SIM_SEED -n $NUM_READS -r $READ_LEN \
  "$FASTA" -o "$TMPDIR/reads.blow5" --paf "$TMPDIR/truth.paf" 2>/dev/null

# 2. Index
"$PANOMAP" index "$GFA" -o "$TMPDIR/drb1.pirx" 2>/dev/null

# 3. Dump node coords
"$PANOMAP" inspect "$TMPDIR/drb1.pirx" --dump-path-coords "$TMPDIR/nodes.tsv" 2>/dev/null

# 4. Map
"$PANOMAP" map --index "$TMPDIR/drb1.pirx" "$TMPDIR/reads.blow5" -o "$TMPDIR/out.gaf" 2>/dev/null

# 5. Eval accuracy
OUTPUT=$("$PYTHON" "$SCRIPTS/eval_canonical.py" "$GFA" "$TMPDIR/nodes.tsv" "$TMPDIR/out.gaf" 2>&1)
echo "$OUTPUT"

# Parse results
CORRECT=$(echo "$OUTPUT" | grep "^Correct:" | grep -oP '\d+')
TOTAL=$(echo "$OUTPUT" | grep "^Total:" | grep -oP '\d+')
WRONG=$(echo "$OUTPUT" | grep "^Wrong:" | grep -oP '\d+')
MAPPED=$(echo "$OUTPUT" | grep "^Mapped:" | grep -oP '\d+')

if [ -z "$TOTAL" ] || [ "$TOTAL" -eq 0 ]; then
  echo "FAIL: no reads evaluated"
  exit 1
fi

RECALL=$((MAPPED * 100 / TOTAL))
MAP_ACC=$((MAPPED > 0 ? CORRECT * 100 / MAPPED : 0))

echo ""
echo "Result: recall=$RECALL% mapping_accuracy=$MAP_ACC% ($CORRECT/$TOTAL correct, $WRONG wrong)"

if [ "$RECALL" -lt "$MIN_ACCURACY" ]; then
  echo "FAIL: recall $RECALL% < ${MIN_ACCURACY}% threshold"
  exit 1
fi

echo "PASS"
