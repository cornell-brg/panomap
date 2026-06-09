#!/usr/bin/env python3
"""
GFA path length calculator.

Parses a GFA file, walks each path, and computes path length by summing
node sequence lengths. Supports both standard GFA (DNA sequences) and
PANOMAP's signal dump GFAs (comma-separated numeric values).

Usage:
    python gfa_path_lengths.py input.gfa
    python gfa_path_lengths.py input.gfa --output lengths.tsv
    python gfa_path_lengths.py input.gfa --verbose  # show per-node details
    python gfa_path_lengths.py *.gfa --summary      # compare multiple GFAs
"""

import argparse
import sys
import re
from pathlib import Path
from typing import Dict, List, Tuple, Optional


def get_segment_length(sequence: str, tags: Dict[str, str]) -> int:
    """
    Get segment length, handling both DNA sequences and signal dumps.

    For signal dumps (comma-separated values), use LN tag if available,
    otherwise count comma-separated elements.
    """
    # Check if this looks like a signal dump (comma-separated numbers)
    if ',' in sequence and sequence.replace(',', '').replace('.', '').replace('-', '').replace(' ', '').isdigit():
        # Signal dump format - use LN tag or count elements
        if 'LN' in tags:
            return int(tags['LN'])
        return len(sequence.split(','))
    else:
        # Standard DNA sequence
        return len(sequence)


def normalize_original_id(oi_value: str) -> str:
    """
    Normalize original ID from oi:Z: tag for matching with path steps.

    Examples:
        '1+' -> '1'  (strip orientation)
        '1-' -> '1'
    """
    return oi_value.rstrip('+-')


def parse_path_step(step: str) -> Tuple[str, str]:
    """
    Parse a path step and return (segment_key, orientation).

    Handles multiple formats:
        '1+' -> ('1', '+')           # standard GFA
        '1_F+' -> ('1', '+')         # PANOMAP transformed format (_F = forward)
        '1_R+' -> ('1', '+')         # PANOMAP reverse format (_R = reverse)
    """
    step = step.strip()
    if not step:
        return ('', '')

    # Extract orientation (last char)
    orientation = step[-1] if step[-1] in '+-' else '+'
    base = step.rstrip('+-')

    # Handle PANOMAP format: X_F or X_R suffix
    if '_F' in base or '_R' in base:
        # '1_F' -> '1', '123_R' -> '123'
        base = re.sub(r'_[FR]$', '', base)

    return (base, orientation)


def parse_gfa(gfa_path: Path) -> Tuple[Dict[str, Tuple[str, int]], List[Tuple[str, List[str], Optional[str]]], str]:
    """
    Parse GFA file and return:
    - segments: dict mapping segment_key -> (sequence, length)
      where segment_key is either the segment ID or normalized original ID (from oi:Z: tag)
    - paths: list of (path_name, [segment_keys], overlap_string)
    - gfa_type: detected type ('dna', 'raw_signal', 'fuzzy_quant', 'aln_quant', 'unknown')
    """
    segments = {}
    original_id_map = {}  # oi -> seg_id mapping for PANOMAP dumps
    paths = []
    gfa_type = 'unknown'

    with open(gfa_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            fields = line.split('\t')
            record_type = fields[0]

            if record_type == 'S':
                # Segment: S <id> <sequence> [tags...]
                seg_id = fields[1]
                sequence = fields[2]

                # Parse tags (e.g., LN:i:40, st:Z:raw_signal, oi:Z:1+)
                tags = {}
                for tag in fields[3:]:
                    if ':' in tag:
                        parts = tag.split(':')
                        if len(parts) >= 3:
                            tag_name = parts[0]
                            tag_value = parts[2]
                            tags[tag_name] = tag_value

                # Detect GFA type from st tag
                if 'st' in tags:
                    st = tags['st']
                    if st == 'raw_signal':
                        gfa_type = 'raw_signal'
                    elif st == 'fuzzy_quant':
                        gfa_type = 'fuzzy_quant'
                    elif st == 'aln_quant':
                        gfa_type = 'aln_quant'
                elif gfa_type == 'unknown' and sequence.replace(',', '').replace('.', '').replace('-', '').replace(' ', '').isdigit() == False:
                    if all(c in 'ACGTNacgtn' for c in sequence):
                        gfa_type = 'dna'

                seg_len = get_segment_length(sequence, tags)

                # Store by segment ID
                segments[seg_id] = (sequence, seg_len)

                # Also store by normalized original ID if oi:Z: tag present
                if 'oi' in tags:
                    norm_oi = normalize_original_id(tags['oi'])
                    original_id_map[norm_oi] = seg_id
                    # Store directly by original ID for easier path lookup
                    segments[norm_oi] = (sequence, seg_len)

            elif record_type == 'P':
                # Path: P <name> <segment_list> <overlap_list>
                path_name = fields[1]
                segment_list = fields[2]
                overlap_list = fields[3] if len(fields) > 3 else '*'

                # Parse segment list (e.g., "1+,2+,4-,5+" or "1_F+,2_F+,...")
                seg_ids = []
                for seg in segment_list.split(','):
                    seg_key, _ = parse_path_step(seg)
                    if seg_key:
                        seg_ids.append(seg_key)

                paths.append((path_name, seg_ids, overlap_list))

    return segments, paths, gfa_type


def compute_path_length(
    segments: Dict[str, Tuple[str, int]],
    seg_ids: List[str],
) -> Tuple[int, List[Tuple[str, int]]]:
    """
    Compute total path length by summing segment lengths.

    Returns:
    - total_length: sum of segment lengths
    - details: list of (seg_id, length) for each segment
    """
    details = []
    total_length = 0

    for seg_id in seg_ids:
        if seg_id not in segments:
            print(f"Warning: segment '{seg_id}' not found in GFA", file=sys.stderr)
            continue

        _, seg_len = segments[seg_id]
        details.append((seg_id, seg_len))
        total_length += seg_len

    return total_length, details


def process_single_gfa(gfa_path: Path, verbose: bool = False) -> List[Tuple[str, int, int]]:
    """Process a single GFA and return list of (path_name, length, num_segments)."""
    segments, paths, gfa_type = parse_gfa(gfa_path)
    print(f"Parsed {gfa_path.name}: {len(segments)} segments, {len(paths)} paths, type={gfa_type}", file=sys.stderr)

    results = []
    for path_name, seg_ids, overlap_str in paths:
        total_length, details = compute_path_length(segments, seg_ids)
        results.append((path_name, total_length, len(seg_ids)))

        if verbose:
            print(f"\n  Path '{path_name}':", file=sys.stderr)
            for seg_id, seg_len in details[:10]:
                print(f"    {seg_id}: {seg_len}", file=sys.stderr)
            if len(details) > 10:
                print(f"    ... ({len(details) - 10} more segments)", file=sys.stderr)

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Parse GFA file(s) and compute path lengths by summing node sequences."
    )
    parser.add_argument("gfa", type=Path, nargs='+', help="Input GFA file(s)")
    parser.add_argument(
        "--output", "-o", type=Path, default=None,
        help="Output TSV file (default: stdout)"
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true",
        help="Show per-segment details for each path"
    )
    parser.add_argument(
        "--summary", "-s", action="store_true",
        help="Show summary comparison across multiple GFAs"
    )
    args = parser.parse_args()

    # Check all files exist
    for gfa_path in args.gfa:
        if not gfa_path.exists():
            print(f"Error: GFA file not found: {gfa_path}", file=sys.stderr)
            sys.exit(1)

    if args.summary and len(args.gfa) > 1:
        # Summary mode: compare path lengths across multiple GFAs
        all_results = {}  # path_name -> {gfa_name: length}
        gfa_names = []

        for gfa_path in args.gfa:
            gfa_name = gfa_path.stem
            gfa_names.append(gfa_name)
            results = process_single_gfa(gfa_path, verbose=False)
            for path_name, length, _ in results:
                if path_name not in all_results:
                    all_results[path_name] = {}
                all_results[path_name][gfa_name] = length

        # Output summary table
        out = open(args.output, "w") if args.output else sys.stdout
        try:
            out.write("path_name\t" + "\t".join(gfa_names) + "\n")
            for path_name in sorted(all_results.keys()):
                row = [path_name]
                for gfa_name in gfa_names:
                    row.append(str(all_results[path_name].get(gfa_name, 'N/A')))
                out.write("\t".join(row) + "\n")
        finally:
            if args.output:
                out.close()
    else:
        # Single GFA or multiple without summary
        out = open(args.output, "w") if args.output else sys.stdout
        try:
            out.write("path_name\tlength\tnum_segments\n")
            for gfa_path in args.gfa:
                results = process_single_gfa(gfa_path, verbose=args.verbose)
                for path_name, length, num_segs in results:
                    out.write(f"{path_name}\t{length}\t{num_segs}\n")
        finally:
            if args.output:
                out.close()


if __name__ == "__main__":
    main()
