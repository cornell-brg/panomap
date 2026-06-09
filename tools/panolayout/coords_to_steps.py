#!/usr/bin/env python3
"""Read piru-style coords.tsv (node_id, canon_start, canon_end, component_id, ...)
plus a GFA, walk the GFA paths, and emit a path_steps.tsv in panolayout format.
This lets us render piru-generated coords with our shared plotter.
"""
import argparse
import csv
from pathlib import Path


def load_coords(tsv_path, divide_node_id_by=1):
    """node_id -> (canon_start, canon_end). Piru dumps node_id = 2*ordinal (even);
    we divide by 2 to get back to ordinal. Set divide_node_id_by=1 for panolayout's
    dense node IDs.
    """
    out = {}
    with open(tsv_path) as f:
        r = csv.DictReader(f, delimiter='\t')
        for row in r:
            nid = int(row['node_id']) // divide_node_id_by
            out[nid] = (float(row['canon_start']), float(row['canon_end']))
    return out


def gfa_iter_paths_with_lens(gfa):
    """Yield (path_name, [(seg_name, orient), ...]) and a global name->(ordinal, length) map."""
    name_to_ord = {}
    name_to_len = {}
    paths = []
    with open(gfa) as f:
        for line in f:
            if line.startswith('S\t'):
                p = line.rstrip().split('\t', 3)
                if len(p) < 3: continue
                name_to_ord[p[1]] = len(name_to_ord)
                name_to_len[p[1]] = len(p[2])
            elif line.startswith('P\t'):
                p = line.rstrip().split('\t')
                if len(p) < 3: continue
                steps = []
                for s in p[2].split(','):
                    if not s: continue
                    steps.append((s[:-1], s[-1]))
                paths.append((p[1], steps))
            elif line.startswith('W\t'):
                p = line.rstrip().split('\t')
                if len(p) < 7: continue
                pname = f"{p[1]}#{p[2]}#{p[3]}"
                walk = p[6]
                steps = []
                i = 0
                while i < len(walk):
                    c = walk[i]
                    if c not in '><': i += 1; continue
                    orient = '+' if c == '>' else '-'
                    j = i + 1
                    while j < len(walk) and walk[j] not in '><': j += 1
                    steps.append((walk[i+1:j], orient))
                    i = j
                paths.append((pname, steps))
    return paths, name_to_ord, name_to_len


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("gfa")
    ap.add_argument("coords_tsv")
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--panomap", action="store_true",
                    help="Panomap style: node_id is 2*ordinal (even). Else assume dense (panolayout).")
    args = ap.parse_args()

    div = 2 if args.panomap else 1
    coords = load_coords(args.coords_tsv, divide_node_id_by=div)
    paths, n2o, n2l = gfa_iter_paths_with_lens(args.gfa)

    with open(args.out, 'w') as out:
        out.write("path_name\tstep_idx\tnode_id\torient\tcum_bp\tcoord_start\tcoord_end\n")
        for pname, steps in paths:
            cum = 0
            for i, (seg, orient) in enumerate(steps):
                ord_ = n2o.get(seg)
                if ord_ is None: continue
                length = n2l.get(seg, 0)
                cs, ce = coords.get(ord_, (0.0, 0.0))
                out.write(f"{pname}\t{i}\t{ord_}\t{orient}\t{cum}\t{cs}\t{ce}\n")
                cum += length
    print(f"-> {args.out}")


if __name__ == "__main__":
    main()
