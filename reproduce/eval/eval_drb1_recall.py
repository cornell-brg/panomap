#!/usr/bin/env python3
"""Recall summary for a DRB1 leave-one-out run.

Walks <runs_dir>/N{N}/chunks_{C}/on_target.gaf and reports, per (N, chunks):
n_reads, n_mapped, recall (= n_mapped / n_reads).

A read counts as mapped if it has a primary alignment (tp:A:P). All reads
are on-target (simulated from the held-out alleles), so recall is the metric.

Usage:
  eval_drb1_recall.py <runs_dir>
"""
import re
import sys
from pathlib import Path


def parse_gaf(p):
    n = n_mapped = 0
    with open(p) as f:
        for line in f:
            if line.startswith("#"):
                continue
            n += 1
            if "tp:A:P" in line.rstrip().split("\t")[11:]:
                n_mapped += 1
    return n, n_mapped


def main():
    if len(sys.argv) < 2:
        print(__doc__, file=sys.stderr)
        sys.exit(1)
    runs = Path(sys.argv[1])
    if not runs.is_dir():
        print(f"ERR: not a directory: {runs}", file=sys.stderr)
        sys.exit(1)

    Ns = sorted(int(re.search(r"N(\d+)$", d.name).group(1))
                for d in runs.glob("N*") if re.search(r"N\d+$", d.name))
    if not Ns:
        print(f"ERR: no N* dirs under {runs}", file=sys.stderr)
        sys.exit(1)

    print(f"{'N':>6}\t{'chunks':>6}\t{'n_reads':>8}\t{'n_mapped':>8}\t{'recall':>7}")
    for N in Ns:
        for cdir in sorted((runs / f"N{N}").glob("chunks_*"),
                           key=lambda p: int(p.name.split("_")[1])):
            gaf = cdir / "on_target.gaf"
            if not gaf.exists():
                continue
            c = int(cdir.name.split("_")[1])
            n, m = parse_gaf(gaf)
            if n == 0:
                continue
            print(f"{N:>6}\t{c:>6}\t{n:>8}\t{m:>8}\t{m/n:>7.4f}")


if __name__ == "__main__":
    main()
