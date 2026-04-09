#!/usr/bin/env python3
"""Evaluate piru mapping accuracy using canonical 1D coordinates.

Reads GFA for path structure + node base lengths.
Reads node coord TSV (from piru inspect --dump-path-coords) for canonical coords.
Reads GAF for piru mapping results (ci:f:, ce:f:, cc:i: tags).
Compares truth canonical interval against piru canonical interval.

Usage:
  python3 eval_canonical.py <graph.gfa> <node_coords.tsv> <output.gaf>
"""

import sys
import re


def load_gfa_paths(gfa_path):
    """Load GFA: node base lengths and path step lists.

    Returns:
      node_lens: dict node_name -> base length
      paths: dict path_name -> [(node_name, orientation), ...]
    """
    node_lens = {}
    paths = {}
    with open(gfa_path) as f:
        for line in f:
            if line.startswith("S\t"):
                parts = line.strip().split("\t")
                name = parts[1]
                seq = parts[2]
                node_lens[name] = len(seq)
            elif line.startswith("P\t"):
                parts = line.strip().split("\t")
                pname = parts[1]
                steps_str = parts[2]
                steps = []
                for s in steps_str.split(","):
                    orient = s[-1]  # + or -
                    nname = s[:-1]
                    steps.append((nname, orient))
                paths[pname] = steps
    return node_lens, paths


def load_node_coords(tsv_path):
    """Load node canonical coords from inspect dump.

    Returns: dict fwd_node_id -> (canon_start, canon_end, component_id)
    """
    coords = {}
    with open(tsv_path) as f:
        next(f)  # header
        for line in f:
            parts = line.strip().split("\t")
            nid = int(parts[0])
            cs = float(parts[1])
            ce = float(parts[2])
            cc = int(parts[3])
            coords[nid] = (cs, ce, cc)
    return coords


def load_node_id_map(gfa_path):
    """Build map from GFA node name to directional graph node IDs.

    simpleExpand creates: forward = orig_idx*2, reverse = orig_idx*2+1.
    Node order = order of S lines in GFA.
    """
    name_to_fwd_id = {}
    idx = 0
    with open(gfa_path) as f:
        for line in f:
            if line.startswith("S\t"):
                name = line.strip().split("\t")[1]
                name_to_fwd_id[name] = idx * 2
                idx += 1
    return name_to_fwd_id


def pos_to_canonical(base_pos, path_steps, node_lens, node_id_map, node_coords):
    """Convert a base-space position on a path to canonical 1D coordinate.

    Walks path steps, accumulating base lengths, until we find the node
    containing base_pos. Then interpolates within that node.

    Returns (canonical_coord, component_id) or None.
    """
    cumulative = 0
    last_result = None
    for nname, orient in path_steps:
        nlen = node_lens.get(nname, 0)
        if nlen == 0:
            continue
        fwd_id = node_id_map.get(nname)
        if fwd_id is None:
            cumulative += nlen
            continue
        cs, ce, cc = node_coords.get(fwd_id, (0, 0, 0))

        if cumulative + nlen >= base_pos:
            # base_pos falls within (or at boundary of) this node
            offset_in_node = base_pos - cumulative
            frac = offset_in_node / nlen if nlen > 0 else 0

            # For reverse orientation, offset is from other end
            if orient == "-":
                frac = 1.0 - frac
            canon = cs + frac * (ce - cs)
            return canon, cc
        cumulative += nlen
        # Track last node in case base_pos is at/past end
        last_result = (ce if orient == "+" else cs, cc)
    # At or past end of path: return last node's end
    return last_result


def overlap_frac(s1, e1, s2, e2):
    """Reciprocal overlap fraction (using min/max to handle inverted intervals)."""
    a_lo, a_hi = min(s1, e1), max(s1, e1)
    b_lo, b_hi = min(s2, e2), max(s2, e2)
    ov = max(0, min(a_hi, b_hi) - max(a_lo, b_lo))
    if ov == 0:
        return 0.0
    la = a_hi - a_lo
    lb = b_hi - b_lo
    if la <= 0 or lb <= 0:
        return 0.0
    return min(ov / la, ov / lb)


def parse_gaf_tags(fields):
    tags = {}
    for f in fields[12:]:
        if ":" in f:
            key = f[:2]
            tags[key] = f
    return tags


def main():
    if len(sys.argv) < 4:
        print(f"Usage: {sys.argv[0]} <graph.gfa> <node_coords.tsv> <output.gaf>")
        sys.exit(1)

    gfa_path = sys.argv[1]
    coords_path = sys.argv[2]
    gaf_path = sys.argv[3]

    node_lens, gfa_paths = load_gfa_paths(gfa_path)
    node_coords = load_node_coords(coords_path)
    node_id_map = load_node_id_map(gfa_path)

    total = 0
    correct = 0
    wrong = 0
    unmapped = 0
    no_truth = 0
    details = []

    with open(gaf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            tags = parse_gaf_tags(fields)
            if "tp" not in tags:
                continue

            is_primary = tags["tp"] == "tp:A:P"
            is_unmapped = tags["tp"] == "tp:A:U"
            if not is_primary and not is_unmapped:
                continue

            read_id = fields[0]
            parts = read_id.split("!")
            if len(parts) < 5:
                continue
            total += 1

            truth_contig = parts[1]
            truth_start = int(parts[2])
            truth_end = int(parts[3])

            if is_unmapped:
                unmapped += 1
                details.append(f"UNMAPPED   {read_id}")
                continue

            if "ci" not in tags or "ce" not in tags or "cc" not in tags:
                no_truth += 1
                continue

            piru_ci = float(tags["ci"].split(":")[2])
            piru_ce = float(tags["ce"].split(":")[2])
            piru_cc = int(tags["cc"].split(":")[2])

            # Get path steps for truth contig (always forward path)
            path_steps = gfa_paths.get(truth_contig)
            if path_steps is None:
                no_truth += 1
                details.append(f"NO_PATH   {read_id}  contig={truth_contig}")
                continue

            result_s = pos_to_canonical(truth_start, path_steps, node_lens, node_id_map, node_coords)
            result_e = pos_to_canonical(truth_end, path_steps, node_lens, node_id_map, node_coords)

            if result_s is None or result_e is None:
                no_truth += 1
                details.append(f"NO_CANON  {read_id}  {truth_contig}:{truth_start}-{truth_end}")
                continue

            truth_ci, truth_cc_s = result_s
            truth_ce, truth_cc_e = result_e
            mapq = int(fields[11])

            # Compare canonical intervals (ignore component for now, check overlap)
            ov = overlap_frac(truth_ci, truth_ce, piru_ci, piru_ce)

            # Truth node walk: nodes covering [truth_start, truth_end] on the path
            truth_nodes = []
            cum = 0
            for nname, orient in path_steps:
                nlen = node_lens.get(nname, 0)
                if cum + nlen > truth_start and cum < truth_end:
                    truth_nodes.append(nname)
                cum += nlen
                if cum >= truth_end:
                    break
            truth_walk = ",".join(truth_nodes[:10])
            if len(truth_nodes) > 10:
                truth_walk += f"...({len(truth_nodes)} nodes)"

            # Piru node walk from GAF col 6
            piru_walk = fields[5][:80]
            if len(fields[5]) > 80:
                piru_walk += "..."

            if ov > 0:
                correct += 1
                details.append(
                    f"OK        {read_id}  "
                    f"truth=[{min(truth_ci,truth_ce):.0f},{max(truth_ci,truth_ce):.0f}]  "
                    f"piru=[{min(piru_ci,piru_ce):.0f},{max(piru_ci,piru_ce):.0f}]  "
                    f"ov={ov:.1%}  cc_t={truth_cc_s}  cc_p={piru_cc}  mq={mapq}"
                )
            else:
                wrong += 1
                details.append(
                    f"WRONG     {read_id}  "
                    f"truth=[{min(truth_ci,truth_ce):.0f},{max(truth_ci,truth_ce):.0f}]  "
                    f"piru=[{min(piru_ci,piru_ce):.0f},{max(piru_ci,piru_ce):.0f}]  "
                    f"ov=0  cc_t={truth_cc_s}  cc_p={piru_cc}  mq={mapq}\n"
                    f"          truth_nodes: {truth_walk}\n"
                    f"          piru_walk:   {piru_walk}"
                )

    pct = correct / total * 100 if total > 0 else 0
    print(f"Total: {total}")
    print(f"Correct: {correct} ({pct:.1f}%)")
    print(f"Wrong: {wrong}")
    print(f"Unmapped: {unmapped}")
    print(f"No truth: {no_truth}")
    print()

    # Print non-OK
    for d in details:
        if not d.startswith("OK"):
            print(d)


if __name__ == "__main__":
    main()
