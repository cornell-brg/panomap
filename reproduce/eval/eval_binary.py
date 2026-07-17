#!/usr/bin/env python3
"""Binary classification metrics from panomap GAF output (covid / yeast).

Reads on-target and off-target panomap GAFs and emits a TSV row with
TP/FP/FN/TN, precision, recall, specificity, F1.

A read is "kept" (positive call) if it has a primary alignment (tp:A:P).
  TP = on-target reads kept        FP = off-target reads kept
  FN = on-target total  - TP       TN = off-target total - FP

  precision   = TP / (TP + FP)
  recall      = TP / (TP + FN)
  specificity = TN / (TN + FP)
  F1          = 2 * precision * recall / (precision + recall)

Usage:
  eval_binary.py --on ON.gaf --off OFF.gaf \
                 --on-total N --off-total N \
                 --name RUN_NAME --chunks N
  eval_binary.py --header        # print TSV header only
"""

import argparse


def count_kept(path):
    """Count distinct reads with a primary alignment (tp:A:P) in a panomap GAF."""
    kept = set()
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 12:
                continue
            if 'tp:A:P' in cols[12:]:
                kept.add(cols[0])
    return len(kept)


HEADER = "name\tchunks\tTP\tFP\tFN\tTN\tprecision\trecall\tspecificity\tF1"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--on')
    ap.add_argument('--off')
    ap.add_argument('--on-total', type=int)
    ap.add_argument('--off-total', type=int)
    ap.add_argument('--name')
    ap.add_argument('--chunks', type=int)
    ap.add_argument('--header', action='store_true', help='Print TSV header only')
    args = ap.parse_args()

    if args.header:
        print(HEADER)
        return

    required = ['on', 'off', 'on_total', 'off_total', 'name', 'chunks']
    missing = [r for r in required if getattr(args, r) is None]
    if missing:
        ap.error(f"missing required args: {missing}")

    tp = count_kept(args.on)
    fp = count_kept(args.off)
    fn = args.on_total - tp
    tn = args.off_total - fp

    prec = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    spec = tn / (tn + fp) if (tn + fp) > 0 else 0.0
    f1 = 2 * prec * recall / (prec + recall) if (prec + recall) > 0 else 0.0

    print(f"{args.name}\t{args.chunks}\t{tp}\t{fp}\t{fn}\t{tn}\t"
          f"{prec:.3f}\t{recall:.3f}\t{spec:.3f}\t{f1:.3f}", flush=True)


if __name__ == '__main__':
    main()
