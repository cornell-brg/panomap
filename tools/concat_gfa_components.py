#!/usr/bin/env python3
"""concat_gfa_components.py -- concatenate per-species GFAs into a
multi-component multi-species graph.

Used to build the multi-species pangenome graphs. Both pggb and
Minigraph-Cactus per-species GFAs are merged the same way.

Each input GFA has its own numeric segment ID space (1..N_segments).
Output GFA assigns each species a non-overlapping ID range by adding a
per-species offset to all S, L, and P line segment references.

Usage:
  concat_gfa_components.py -o OUT.gfa SPECIES1=PATH1 SPECIES2=PATH2 ...

Example:
  concat_gfa_components.py -o zymo-7sp-8.gfa \\
      ecoli=ecoli-8.gfa \\
      bsubtilis=bsubtilis-8.gfa \\
      ...

Notes:
- Path names (col 2 of P lines) are kept as-is (NCBI accessions are
  globally unique).
- Walk (W) lines are translated similarly if present.
- Header (H) line emitted once at the top from the first input.
"""

import argparse
import sys
from pathlib import Path


def rewrite_path_segments(seg_csv, offset):
    """Rewrite a P-line seg list 'id1+,id2-,id3+,...' with offset added."""
    parts = seg_csv.split(",")
    out = []
    for p in parts:
        if not p:
            continue
        orient = p[-1]
        sid = p[:-1]
        out.append(f"{int(sid) + offset}{orient}")
    return ",".join(out)


def rewrite_walk_segments(walk_str, offset):
    """Rewrite a W-line walk '>id1<id2>id3...' with offset added."""
    out = []
    i = 0
    while i < len(walk_str):
        marker = walk_str[i]
        if marker not in ("<", ">"):
            raise ValueError(f"Bad walk char {marker!r} at pos {i}")
        i += 1
        start = i
        while i < len(walk_str) and walk_str[i] not in ("<", ">"):
            i += 1
        sid = walk_str[start:i]
        out.append(f"{marker}{int(sid) + offset}")
    return "".join(out)


def process_species(in_path, offset, species_tag, header_done, out_fh):
    """Stream-rewrite one input GFA. Returns max numeric segment id seen.

    Handles two GFA dialects:
    - Numeric seg IDs (typical for N>=2): apply offset.
    - String seg IDs (N=1, or anywhere accession-named): prefix with
      species_tag for cross-species disambiguation.
    """
    max_seg = 0
    # Probe first S line to decide dialect
    string_mode = False
    with open(in_path) as f:
        for probe in f:
            if probe.startswith("S\t"):
                first_id = probe.split("\t", 2)[1]
                try:
                    int(first_id)
                except ValueError:
                    string_mode = True
                break

    def xform_seg(sid):
        if string_mode:
            return f"{species_tag}_{sid}"
        return str(int(sid) + offset)

    with open(in_path) as f:
        for line in f:
            if not line:
                continue
            tag = line[0]
            if tag == "H":
                if not header_done[0]:
                    out_fh.write(line)
                    header_done[0] = True
                continue
            cols = line.rstrip("\n").split("\t")
            if tag == "S":
                if not string_mode:
                    max_seg = max(max_seg, int(cols[1]) + offset)
                cols[1] = xform_seg(cols[1])
                out_fh.write("\t".join(cols) + "\n")
            elif tag == "L":
                cols[1] = xform_seg(cols[1])
                cols[3] = xform_seg(cols[3])
                out_fh.write("\t".join(cols) + "\n")
            elif tag == "P":
                cols[1] = f"{species_tag}#{cols[1]}"
                if string_mode:
                    # Walk format: id1+,id2-,... but ids are strings
                    parts = cols[2].split(",")
                    new_parts = []
                    for p in parts:
                        if not p:
                            continue
                        orient = p[-1]
                        sid = p[:-1]
                        new_parts.append(f"{xform_seg(sid)}{orient}")
                    cols[2] = ",".join(new_parts)
                else:
                    cols[2] = rewrite_path_segments(cols[2], offset)
                out_fh.write("\t".join(cols) + "\n")
            elif tag == "W":
                cols[1] = f"{species_tag}#{cols[1]}"
                if string_mode:
                    # Walk: <id1>id2<id3... -- strings
                    out_walk = []
                    i = 0
                    s = cols[6]
                    while i < len(s):
                        marker = s[i]
                        if marker not in ("<", ">"):
                            raise ValueError(f"Bad walk char {marker!r}")
                        i += 1
                        start = i
                        while i < len(s) and s[i] not in ("<", ">"):
                            i += 1
                        sid = s[start:i]
                        out_walk.append(f"{marker}{xform_seg(sid)}")
                    cols[6] = "".join(out_walk)
                else:
                    cols[6] = rewrite_walk_segments(cols[6], offset)
                out_fh.write("\t".join(cols) + "\n")
            else:
                out_fh.write(line)
    return max_seg


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--output", required=True, help="Output GFA path")
    ap.add_argument("species_paths", nargs="+",
                    help="Per-species inputs as SPECIES=PATH")
    args = ap.parse_args()

    inputs = []
    for sp in args.species_paths:
        if "=" not in sp:
            sys.exit(f"error: bad input '{sp}', expected SPECIES=PATH")
        species, path = sp.split("=", 1)
        if not Path(path).is_file():
            sys.exit(f"error: missing {path}")
        inputs.append((species, path))

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    cumulative_offset = 0
    header_done = [False]

    print(f"Concatenating {len(inputs)} species into {out_path}")
    with out_path.open("w") as out_fh:
        for species, in_path in inputs:
            print(f"  + {species:18s} (offset={cumulative_offset:>10d}) <- {in_path}")
            max_seg = process_species(
                in_path, cumulative_offset, species, header_done, out_fh)
            cumulative_offset = max_seg + 1
    print(f"Done. Final cumulative_offset = {cumulative_offset:,}")
    print(f"Output: {out_path} ({out_path.stat().st_size / 1024**2:.1f} MB)")


if __name__ == "__main__":
    main()
