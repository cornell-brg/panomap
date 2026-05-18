#!/bin/bash
# Index-format oracle test (dev-112).
#
# Builds known graphs at the current piru version, maps a fixed read set
# with -t 1, sorts the GAF by read_id, and diffs against pre-committed
# fixtures in tests/fixtures/oracle_v2/. Any byte-level divergence fails.
#
# Purpose: serialization-only changes (drop edges / drop path steps /
# isReverse derivation) must not perturb chain output. This is the
# regression gate for dev-112-index-slim-v3.
#
# Usage:
#   test_index_oracle.sh <piru_binary> <workspace_root>           # diff mode
#   test_index_oracle.sh <piru_binary> <workspace_root> --gen     # regenerate fixtures
#
# Fixtures: tests/fixtures/oracle_v2/<dataset>-N<N>.gaf
# Inputs reside in workspace (data/eval/...). Test skips a tier if inputs
# are missing on the current host (e.g. clean clone with no smoke data).

set -euo pipefail

PIRU="${1:?Usage: $0 <piru_binary> <workspace_root> [--gen]}"
WORKSPACE="${2:?Usage: $0 <piru_binary> <workspace_root> [--gen]}"
MODE="${3:-test}"  # 'test' (default) or '--gen'

FIXTURE_DIR="$(cd "$(dirname "$0")/.." && pwd)/fixtures/oracle_v2"
mkdir -p "$FIXTURE_DIR"

TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

COVID_GFA_DIR="$WORKSPACE/data/eval/sars-cov-2/chronological/pggb"
COVID_READS="$WORKSPACE/experiments/eval/build-rawhash2-smoke/reads/covid-interartic-r9-40k_split/covid-interartic-r9-40k_0.blow5"
YEAST_GFA_DIR="$WORKSPACE/data/eval/yeast"
YEAST_READS="$WORKSPACE/experiments/eval/build-yeast-piru-signal-landmark/smoke_reads/yeast-r9-40k_split/yeast-r9-40k_0.blow5"

EXIT=0

run_one() {
  local ds=$1 N=$2 gfa=$3 reads=$4
  local fixture="$FIXTURE_DIR/${ds}-N${N}.gaf"
  local idx="$TMPDIR/${ds}-N${N}.pirx"
  local gaf="$TMPDIR/${ds}-N${N}.gaf"

  if [ ! -f "$gfa" ] || [ ! -f "$reads" ]; then
    echo "  SKIP ${ds}-N${N}: missing inputs"
    return 0
  fi

  echo "  ${ds}-N${N}: build index"
  "$PIRU" index "$gfa" --model r9.4 -o "$idx" 2>/dev/null

  echo "  ${ds}-N${N}: map -t 1"
  "$PIRU" map --index "$idx" --chainer path-chain -t 1 "$reads" -o "$gaf.raw" 2>/dev/null

  # Sort GAF (excluding header) by read_id for stable comparison.
  # Strip the dt:f: tag -- per-read processing time is non-deterministic
  # (microsecond wall-clock variance between runs). All other tags are
  # deterministic functions of the input and the chain math, which is
  # what the oracle is gating on.
  head -1 "$gaf.raw" > "$gaf"
  tail -n +2 "$gaf.raw" | sed -E 's/\tdt:f:[^\t]*//g' | sort -k1,1 >> "$gaf"

  if [ "$MODE" = "--gen" ]; then
    cp "$gaf" "$fixture"
    echo "    -> generated $fixture ($(wc -l < $fixture) lines)"
  else
    if [ ! -f "$fixture" ]; then
      echo "    MISSING fixture $fixture (run with --gen first)"
      EXIT=1
      return 0
    fi
    if diff -q "$gaf" "$fixture" >/dev/null 2>&1; then
      echo "    OK ($(wc -l < $gaf) lines match)"
    else
      echo "    FAIL: diff vs $fixture"
      diff "$gaf" "$fixture" | head -20
      EXIT=1
    fi
  fi
}

echo "=== Tier 1: covid (N=1,2,4,8) ==="
for N in 1 2 4 8; do
  run_one covid "$N" "$COVID_GFA_DIR/covid-pggb-${N}.gfa" "$COVID_READS"
done

echo
echo "=== Tier 2: yeast (N=1,4) ==="
for N in 1 4; do
  run_one yeast "$N" "$YEAST_GFA_DIR/yeast-pggb-${N}.gfa" "$YEAST_READS"
done

if [ $EXIT -ne 0 ]; then
  echo
  echo "ORACLE FAILED -- chain output diverges from v2 baseline."
  echo "If this is intentional, regenerate with: $0 $PIRU $WORKSPACE --gen"
  exit $EXIT
fi
echo
echo "ORACLE PASSED"
