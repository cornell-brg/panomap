#!/bin/bash
# Integration test: COVID linear reference with RH2 viral preset params
# Maps 1000 real r9.4 reads, compares against pre-generated minimap2 truth PAF.

set -euo pipefail

PIRU="${1:?Usage: $0 <piru_binary> <workspace_root>}"
WORKSPACE="${2:?Usage: $0 <piru_binary> <workspace_root>}"

DATA="$WORKSPACE/data/benchmark/covid-r94"
GFA="$DATA/covid.gfa"
READS="$DATA/covid_reads_1k.blow5"
TRUTH_PAF="$DATA/covid_1k_mm2.paf"

MIN_RECALL=85
MIN_MAP_ACC=85

for f in "$GFA" "$READS" "$TRUTH_PAF"; do
  if [ ! -f "$f" ]; then echo "SKIP: missing $f"; exit 0; fi
done

TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Index with r9.4, seed-k 6 (viral preset)
"$PIRU" index "$GFA" --model r9.4 --seed-k 6 -o "$TMPDIR/covid.pirx" 2>/dev/null

# Map with viral chaining params
"$PIRU" map --index "$TMPDIR/covid.pirx" "$READS" \
  --chain-bw 100 --chain-max-dist 500 --chain-min-score 10 \
  --chain-pen-gap 1.2 --chain-pen-skip 0.3 --max-chunks 5 \
  -o "$TMPDIR/out.gaf" 2>/dev/null

python3 -c "
import sys

mm2 = {}
with open('$TRUTH_PAF') as f:
    for line in f:
        p = line.strip().split('\t')
        if p[0] not in mm2:
            mm2[p[0]] = (int(p[7]), int(p[8]))

correct = wrong = unmapped = 0
with open('$TMPDIR/out.gaf') as f:
    for line in f:
        if line.startswith('#'): continue
        p = line.strip().split('\t')
        rid = p[0]
        if rid not in mm2: continue
        tags = {t[:2]: t for t in p[12:] if ':' in t}
        if 'tp' not in tags: continue
        if tags['tp'] == 'tp:A:U':
            unmapped += 1
            continue
        if tags['tp'] != 'tp:A:P': continue
        if p[7] == '*':
            unmapped += 1
            continue
        ps, pe = int(p[7]), int(p[8])
        ms, me = mm2[rid]
        ov = max(0, min(pe, me) - max(ps, ms))
        if ov > 0: correct += 1
        else: wrong += 1

total = correct + wrong + unmapped
mapped = correct + wrong
recall = mapped * 100 // total if total > 0 else 0
map_acc = correct * 100 // mapped if mapped > 0 else 0
print(f'Total: {total}')
print(f'Mapped: {mapped}')
print(f'Correct: {correct}')
print(f'Wrong: {wrong}')
print(f'Unmapped: {unmapped}')
print(f'Recall: {recall}%')
print(f'Mapping_accuracy: {map_acc}%')
" | tee "$TMPDIR/eval.txt"

RECALL=$(grep "^Recall:" "$TMPDIR/eval.txt" | grep -oP '\d+')
MAP_ACC=$(grep "^Mapping_accuracy:" "$TMPDIR/eval.txt" | grep -oP '\d+')
CORRECT=$(grep "^Correct:" "$TMPDIR/eval.txt" | grep -oP '\d+')
TOTAL=$(grep "^Total:" "$TMPDIR/eval.txt" | grep -oP '\d+')
WRONG=$(grep "^Wrong:" "$TMPDIR/eval.txt" | grep -oP '\d+')

echo ""
echo "Result: recall=$RECALL% mapping_accuracy=$MAP_ACC% ($CORRECT/$TOTAL correct, $WRONG wrong) [covid r9.4 rh2-k6 viral]"

FAIL=0
if [ "$RECALL" -lt "$MIN_RECALL" ]; then
  echo "FAIL: recall $RECALL% < ${MIN_RECALL}%"
  FAIL=1
fi
if [ "$MAP_ACC" -lt "$MIN_MAP_ACC" ]; then
  echo "FAIL: mapping_accuracy $MAP_ACC% < ${MIN_MAP_ACC}%"
  FAIL=1
fi

if [ "$FAIL" -eq 1 ]; then exit 1; fi
echo "PASS"
