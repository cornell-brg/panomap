#!/usr/bin/env python3
"""Convert a FASTA file to a trivial GFA1 graph.

Each sequence becomes a single node (S line) with a path (P line) containing
just that node. This is the minimal GFA that panomap can index for linear
reference benchmarking.

Usage:
    fa2gfa.py input.fa > output.gfa
    fa2gfa.py input.fa -o output.gfa
"""

import argparse
import sys


def fasta_iter(fh):
    """Yield (name, seq) tuples from a FASTA file handle."""
    name, parts = None, []
    for line in fh:
        line = line.rstrip("\n")
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(parts)
            name = line[1:].split()[0]
            parts = []
        else:
            parts.append(line)
    if name is not None:
        yield name, "".join(parts)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("fasta", help="Input FASTA file")
    parser.add_argument("-o", "--output", default="-",
        help="Output GFA file (default: stdout)")
    args = parser.parse_args()

    out = sys.stdout if args.output == "-" else open(args.output, "w")

    out.write("H\tVN:Z:1.0\n")
    for name, seq in fasta_iter(open(args.fasta)):
        out.write(f"S\t{name}\t{seq}\n")
        out.write(f"P\t{name}\t{name}+\t*\n")

    if out is not sys.stdout:
        out.close()


if __name__ == "__main__":
    main()
