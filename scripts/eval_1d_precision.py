#!/usr/bin/env python3
"""
Evaluate chain location accuracy using 1D canonical coordinates.

For each mapped read:
  1. Parse read ID to get source haplotype + interval (simulation ground truth)
  2. Look up source interval's 1D range via pre-computed .pira annotation
  3. Compare chain's 1D ref range against expected 1D range
  4. Report overlap statistics

Requires: annotated .pira v2 file matching the reads' source intervals.

Usage:
  # Step 1: Generate BED from read IDs
  python3 eval_1d_precision.py extract-bed --gaf on_target.gaf -o reads.bed

  # Step 2: Annotate with panomap
  panomap annotate --bed reads.bed --1d-coords-file coords.tsv --mode strict graph.gfa -o reads.pira

  # Step 3: Evaluate
  python3 eval_1d_precision.py eval --gaf on_target.gaf --pira reads.pira
"""

import argparse
import sys
from collections import defaultdict
from dataclasses import dataclass


@dataclass
class PiraRegion:
    label: str
    start: float
    end: float


def parse_read_source(read_id: str):
    """Extract source path and interval from read ID.
    Format: S1_N!path_name!start!end!strand
    Returns (path_name, start, end) or None if unparseable.
    """
    parts = read_id.split("!")
    if len(parts) < 4:
        return None
    path_name = parts[1]
    try:
        start = int(parts[2])
        end = int(parts[3])
    except ValueError:
        return None
    return (path_name, start, end)


def cmd_extract_bed(args):
    """Extract BED entries from GAF read IDs."""
    seen = set()
    entries = []

    with open(args.gaf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            read_id = line.split("\t")[0]
            if read_id in seen:
                continue
            seen.add(read_id)

            source = parse_read_source(read_id)
            if source is None:
                continue
            path_name, start, end = source
            entries.append((path_name, start, end, read_id))

    with open(args.output, "w") as f:
        for path_name, start, end, read_id in entries:
            f.write(f"{path_name}\t{start}\t{end}\t{read_id}\n")

    print(f"Extracted {len(entries)} BED entries to {args.output}")


def parse_pira(path: str) -> dict[str, PiraRegion]:
    """Parse .pira v2, return dict keyed by region label."""
    regions = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            label = parts[0]
            if parts[1] == "*" or parts[2] == "*":
                continue
            regions[label] = PiraRegion(
                label=label,
                start=float(parts[1]),
                end=float(parts[2]),
            )
    return regions


def intervals_overlap(a_start, a_end, b_start, b_end):
    return a_start < b_end and b_start < a_end


def overlap_fraction(a_start, a_end, b_start, b_end):
    a_len = a_end - a_start
    if a_len <= 0:
        return 0.0
    ov_start = max(a_start, b_start)
    ov_end = min(a_end, b_end)
    return max(0, ov_end - ov_start) / a_len


def cmd_eval(args):
    """Evaluate chain location accuracy."""
    regions = parse_pira(args.pira)
    if not regions:
        print("ERROR: no regions in .pira", file=sys.stderr)
        return 1

    # Parse GAF, match each read to its expected region
    correct = 0
    incorrect = 0
    no_region = 0
    total = 0
    correct_scores = []
    incorrect_scores = []

    with open(args.gaf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 12:
                continue

            # Only primary
            if args.primary_only and "tp:A:P" not in line:
                continue

            read_id = cols[0]
            chain_start = int(cols[7])
            chain_end = int(cols[8])

            # Get chain score
            chain_score = 0
            for col in cols[12:]:
                if col.startswith("cs:i:"):
                    chain_score = int(col[5:])

            # Find this read's expected region
            source = parse_read_source(read_id)
            if source is None:
                no_region += 1
                continue

            path_name, src_start, src_end = source
            # Build the label as annotate would: "path:start-end"
            region_label = f"{path_name}:{src_start}-{src_end}"

            if region_label not in regions:
                no_region += 1
                continue

            region = regions[region_label]
            total += 1

            if intervals_overlap(chain_start, chain_end, region.start, region.end):
                correct += 1
                correct_scores.append(chain_score)
            else:
                incorrect += 1
                incorrect_scores.append(chain_score)

    accuracy = correct / total if total > 0 else 0.0

    print(f"Chain location accuracy")
    print(f"  Reads evaluated: {total}")
    print(f"  Correct (chain overlaps expected 1D interval): {correct}")
    print(f"  Incorrect (chain outside expected interval):   {incorrect}")
    print(f"  No region found (skipped):                     {no_region}")
    print(f"  Accuracy: {accuracy:.3f}")
    print()
    if correct_scores:
        print(f"  Correct chain scores:   min={min(correct_scores)} median={sorted(correct_scores)[len(correct_scores)//2]} max={max(correct_scores)} mean={sum(correct_scores)/len(correct_scores):.0f}")
    if incorrect_scores:
        print(f"  Incorrect chain scores: min={min(incorrect_scores)} median={sorted(incorrect_scores)[len(incorrect_scores)//2]} max={max(incorrect_scores)} mean={sum(incorrect_scores)/len(incorrect_scores):.0f}")

    return 0


def main():
    parser = argparse.ArgumentParser(description="Evaluate chain 1D location accuracy")
    sub = parser.add_subparsers(dest="command")

    # extract-bed
    p_bed = sub.add_parser("extract-bed", help="Extract BED from GAF read IDs")
    p_bed.add_argument("--gaf", required=True)
    p_bed.add_argument("-o", "--output", required=True)

    # eval
    p_eval = sub.add_parser("eval", help="Evaluate chain location accuracy")
    p_eval.add_argument("--gaf", required=True)
    p_eval.add_argument("--pira", required=True)
    p_eval.add_argument("--primary-only", action="store_true", default=True)

    args = parser.parse_args()
    if args.command == "extract-bed":
        return cmd_extract_bed(args)
    elif args.command == "eval":
        return cmd_eval(args)
    else:
        parser.print_help()
        return 1


if __name__ == "__main__":
    sys.exit(main())
