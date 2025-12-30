#!/usr/bin/env python3
"""
FASTA sequence length statistics.

Parses a FASTA file and outputs per-sequence length statistics.

Usage:
    python fasta_stats.py input.fa
    python fasta_stats.py input.fa --output stats.tsv
"""

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Tuple


def parse_fasta(fasta_path: Path) -> List[Tuple[str, int]]:
    """
    Parse FASTA file and return list of (sequence_name, length) tuples.
    """
    sequences = []
    current_name = None
    current_length = 0

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                # Save previous sequence if exists
                if current_name is not None:
                    sequences.append((current_name, current_length))
                # Start new sequence
                current_name = line[1:].split()[0]  # Take first word after >
                current_length = 0
            else:
                current_length += len(line)

    # Don't forget the last sequence
    if current_name is not None:
        sequences.append((current_name, current_length))

    return sequences


def main():
    parser = argparse.ArgumentParser(
        description="Parse FASTA file and output per-sequence length statistics."
    )
    parser.add_argument("fasta", type=Path, help="Input FASTA file")
    parser.add_argument(
        "--output", "-o", type=Path, default=None,
        help="Output TSV file (default: stdout)"
    )
    args = parser.parse_args()

    if not args.fasta.exists():
        print(f"Error: FASTA file not found: {args.fasta}", file=sys.stderr)
        sys.exit(1)

    sequences = parse_fasta(args.fasta)

    # Output
    out = open(args.output, "w") if args.output else sys.stdout
    try:
        out.write("sequence_name\tlength\n")
        total_length = 0
        for name, length in sequences:
            out.write(f"{name}\t{length}\n")
            total_length += length

        # Summary to stderr
        print(f"Parsed {len(sequences)} sequences, total {total_length} bp", file=sys.stderr)
    finally:
        if args.output:
            out.close()


if __name__ == "__main__":
    main()
