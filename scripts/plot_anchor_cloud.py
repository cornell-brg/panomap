#!/usr/bin/env python3
"""Plot anchor clouds from panomap trace dumps.

Visualize the anchor cloud (query_pos vs ref_coord) for one or more reads,
highlighting the best chain. Useful for diagnosing why reads map or don't
map -- comparing chain geometry between true matches (colinear diagonal)
and noise (random scatter).

Prerequisites:
    pip install matplotlib numpy

Enable tracing:
    cmake -S repo -B repo/build -DCMAKE_CXX_FLAGS="-DPIRU_TRACE_ENABLED"
    cmake --build repo/build -j8

Generate trace dumps:
    PANOMAP_TRACE_STAGES=0x40 PANOMAP_TRACE_DIR=/tmp/trace \
        panomap map --index ref.pirx reads.blow5 --no-map-filter -o out.gaf

    Stage 0x40 = kChains. Produces per-read files:
        /tmp/trace/6_anchors_<read_id>_chunk<N>
        /tmp/trace/7_decision_<read_id>

Usage:
    # Plot a single trace file:
    python plot_anchor_cloud.py /tmp/trace/6_anchors_S1_1_chunk4

    # Plot all reads in a trace directory (one PNG per read):
    python plot_anchor_cloud.py /tmp/trace/ --all

    # Plot specific reads side by side (comparison mode):
    python plot_anchor_cloud.py /tmp/trace/ \
        --reads S1_1 S1_2 S1_3 --label "human noise" \
        --reads2 M1_1 M1_2 M1_3 --label2 "microbe noise"

    # Use the last chunk (most anchors) by default, or specify:
    python plot_anchor_cloud.py /tmp/trace/ --reads S1_1 --chunk 3
"""

import argparse
import re
import sys
from pathlib import Path
from collections import defaultdict

import numpy as np


def parse_anchor_file(path):
    """Parse a trace 6_anchors file.

    Returns:
        anchors: list of (query_pos, ref_coord, node_id, offset, span)
        chains: list of dicts with 'score', 'anchors' (list of tuples),
                and 'chain_id'
    """
    anchors = []
    chains = []
    current_chain = None

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("# Anchor cloud:"):
                continue
            if line.startswith("# Chains:"):
                continue

            m = re.match(r"# chain (\d+) score=([\d.]+) anchors=(\d+)", line)
            if m:
                current_chain = {
                    "chain_id": int(m.group(1)),
                    "score": float(m.group(2)),
                    "expected_anchors": int(m.group(3)),
                    "anchors": [],
                }
                chains.append(current_chain)
                continue

            parts = line.split("\t")
            if len(parts) >= 5:
                # Check if this is a chain anchor (has CHAIN=N suffix)
                chain_tag = None
                for p in parts[5:]:
                    cm = re.match(r"CHAIN=(\d+)", p)
                    if cm:
                        chain_tag = int(cm.group(1))

                row = (
                    float(parts[0]),  # query_pos
                    float(parts[1]),  # ref_coord
                    int(parts[2]),    # node_id
                    int(parts[3]),    # offset
                    int(parts[4]),    # span / length
                )

                if chain_tag is not None and current_chain is not None:
                    current_chain["anchors"].append(row)
                else:
                    anchors.append(row)

    return anchors, chains


def parse_decision_file(trace_dir, read_id):
    """Find and parse a trace 7_decision file. Returns dict of decision scores."""
    # Try multiple naming patterns
    candidates = [
        trace_dir / f"7_decision_{read_id}",
        trace_dir / f"7_decision_{read_id}.txt",
    ]
    for path in candidates:
        if path.exists():
            with open(path) as f:
                line = f.readline().strip()
            result = {}
            for part in line.split("\t"):
                if "=" in part:
                    k, v = part.split("=", 1)
                    try:
                        result[k] = float(v)
                    except ValueError:
                        result[k] = v
            return result
    return None


def find_trace_files(trace_dir, read_substr=None):
    """Find all 6_anchors files, grouped by read_id.

    Returns: dict of read_id -> list of (chunk_idx, path) sorted by chunk.
    """
    trace_dir = Path(trace_dir)
    files = defaultdict(list)

    for f in sorted(trace_dir.glob("6_anchors_*")):
        name = f.stem if f.suffix else f.name
        # Format: 6_anchors_<read_id>_chunk<N> or 6_anchors_<read_id>_c<N>
        m = re.match(r"6_anchors_(.+?)_(?:chunk|c)(\d+)$", name)
        if not m:
            continue
        rid = m.group(1)
        chunk = int(m.group(2))

        if read_substr and not any(s in rid for s in read_substr):
            continue

        files[rid].append((chunk, f))

    # Sort by chunk index
    for rid in files:
        files[rid].sort(key=lambda x: x[0])

    return dict(files)


def plot_single(ax, anchors, chains, title="", color_chain="green"):
    """Plot one anchor cloud on a matplotlib axes."""
    if not anchors and not chains:
        ax.set_title(title + "\n(no anchors)", fontsize=9)
        return

    # All cloud anchors (gray)
    if anchors:
        qp = [a[0] for a in anchors]
        rc = [a[1] for a in anchors]
        ax.scatter(qp, rc, s=2, alpha=0.1, c="gray",
                   label=f"cloud ({len(anchors)})")

    # Best chain (colored)
    if chains:
        best = chains[0]
        if best["anchors"]:
            bq = [a[0] for a in best["anchors"]]
            br = [a[1] for a in best["anchors"]]
            ax.scatter(bq, br, s=18, alpha=0.8, c=color_chain, zorder=5,
                       label=f"chain (n={len(best['anchors'])}, "
                             f"s={best['score']:.0f})")

            # Connect chain anchors
            order = sorted(range(len(bq)), key=lambda i: bq[i])
            ax.plot([bq[i] for i in order], [br[i] for i in order],
                    c=color_chain, alpha=0.4, linewidth=1, zorder=4)

            # Zoom to chain region only (not full cloud which spans genome)
            q_range = max(bq) - min(bq)
            r_range = max(br) - min(br)
            margin = max(q_range, r_range) * 0.2 + 50
            ax.set_xlim(min(bq) - margin, max(bq) + margin)
            ax.set_ylim(min(br) - margin, max(br) + margin)

    # Title
    score_str = ""
    n_anchors = 0
    if chains and chains[0]["anchors"]:
        score_str = f"score={chains[0]['score']:.0f}"
        n_anchors = len(chains[0]["anchors"])
    ax.set_title(f"{title}\n{score_str}  anchors={n_anchors}",
                 fontsize=9)
    ax.set_xlabel("query_pos", fontsize=8)
    ax.set_ylabel("ref_coord", fontsize=8)
    ax.legend(fontsize=6, loc="upper left")
    ax.tick_params(labelsize=7)


def plot_reads(trace_dir, read_ids, chunk=None, out=None, dpi=150,
               title_prefix="", color="green", label=""):
    """Plot anchor clouds for a list of reads, one panel per read."""
    import matplotlib.pyplot as plt

    trace_dir = Path(trace_dir)
    all_files = find_trace_files(trace_dir, read_ids)

    if not all_files:
        print(f"No trace files found for reads: {read_ids}", file=sys.stderr)
        return

    # Filter to requested reads (preserve order)
    ordered = []
    for rid_substr in read_ids:
        for rid, chunks in all_files.items():
            if rid_substr in rid and rid not in [r for r, _ in ordered]:
                ordered.append((rid, chunks))

    if not ordered:
        ordered = list(all_files.items())

    n = len(ordered)
    fig, axes = plt.subplots(1, n, figsize=(5 * n, 4.5), squeeze=False)
    axes = axes[0]

    for idx, (rid, chunk_files) in enumerate(ordered):
        # Pick chunk: last by default (most anchors), or specified
        if chunk is not None:
            selected = [f for c, f in chunk_files if c == chunk]
            if not selected:
                selected = [chunk_files[-1][1]]
            else:
                selected = [selected[0]]
        else:
            selected = [chunk_files[-1][1]]

        anchors, chains = parse_anchor_file(selected[0])

        # Try to load decision
        decision = parse_decision_file(trace_dir, rid)
        ws = f"ws={decision['w']:.3f}" if decision else ""

        # Short read ID for title
        parts = rid.split("!")
        if len(parts) >= 3:
            short = f"{parts[0]}  {parts[1]}:{parts[2]}"
        else:
            short = rid[:40]

        chunk_num = chunk_files[-1][0] if chunk is None else chunk
        plot_single(axes[idx], anchors, chains,
                    title=f"{short}\nchunk={chunk_num} {ws}",
                    color_chain=color)

    suptitle = label or title_prefix or "Anchor clouds"
    fig.suptitle(suptitle, fontsize=12, y=1.02)
    plt.tight_layout()

    if out:
        fig.savefig(out, dpi=dpi, bbox_inches="tight")
        print(f"Saved: {out}")
    else:
        plt.show()
    plt.close()


def plot_comparison(trace_dir, reads_a, reads_b, label_a="Group A",
                    label_b="Group B", chunk=None, out=None, dpi=150):
    """Plot two groups of reads side by side for comparison.

    Two rows: top row = group A, bottom row = group B.
    """
    import matplotlib.pyplot as plt

    trace_dir = Path(trace_dir)
    all_files = find_trace_files(trace_dir)

    def resolve(read_substrs):
        result = []
        for sub in read_substrs:
            for rid, chunks in all_files.items():
                if sub in rid and rid not in [r for r, _ in result]:
                    result.append((rid, chunks))
                    break
        return result

    group_a = resolve(reads_a)
    group_b = resolve(reads_b)

    n_cols = max(len(group_a), len(group_b))
    if n_cols == 0:
        print("No matching reads found", file=sys.stderr)
        return

    fig, axes = plt.subplots(2, n_cols, figsize=(5 * n_cols, 9),
                             squeeze=False)

    for row, group, label, color in [
        (0, group_a, label_a, "green"),
        (1, group_b, label_b, "red"),
    ]:
        for col, (rid, chunk_files) in enumerate(group):
            if chunk is not None:
                sel = [f for c, f in chunk_files if c == chunk]
                fpath = sel[0] if sel else chunk_files[-1][1]
            else:
                fpath = chunk_files[-1][1]

            anchors, chains = parse_anchor_file(fpath)
            decision = parse_decision_file(trace_dir, rid)
            ws = f"ws={decision['w']:.3f}" if decision else ""

            parts = rid.split("!")
            short = f"{parts[0]} {parts[1]}:{parts[2]}" if len(parts) >= 3 else rid[:40]
            chunk_num = chunk_files[-1][0] if chunk is None else chunk

            plot_single(axes[row][col], anchors, chains,
                        title=f"{short}\nchunk={chunk_num} {ws}",
                        color_chain=color)

        # Label remaining empty panels
        for col in range(len(group), n_cols):
            axes[row][col].set_visible(False)

        # Row label
        axes[row][0].set_ylabel(f"{label}\nref_coord", fontsize=9)

    fig.suptitle(f"{label_a} vs {label_b}", fontsize=13, y=1.02)
    plt.tight_layout()

    if out:
        fig.savefig(out, dpi=dpi, bbox_inches="tight")
        print(f"Saved: {out}")
    else:
        plt.show()
    plt.close()


def main():
    p = argparse.ArgumentParser(
        description="Plot anchor clouds from panomap trace dumps",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("trace_path",
                   help="Trace directory or single 6_anchors file")
    p.add_argument("--reads", nargs="+",
                   help="Read ID substrings to plot (group A)")
    p.add_argument("--reads2", nargs="+",
                   help="Read ID substrings for comparison (group B)")
    p.add_argument("--label", default="",
                   help="Label for group A (default: auto)")
    p.add_argument("--label2", default="",
                   help="Label for group B")
    p.add_argument("--all", action="store_true",
                   help="Plot all reads (one PNG per read)")
    p.add_argument("--chunk", type=int, default=None,
                   help="Chunk index to plot (default: last/highest)")
    p.add_argument("--out", help="Output path (PNG/PDF/SVG)")
    p.add_argument("--dpi", type=int, default=150)
    args = p.parse_args()

    trace_path = Path(args.trace_path)

    # Single file mode
    if trace_path.is_file():
        import matplotlib.pyplot as plt
        anchors, chains = parse_anchor_file(trace_path)
        fig, ax = plt.subplots(1, 1, figsize=(6, 5))
        plot_single(ax, anchors, chains, title=trace_path.stem)
        plt.tight_layout()
        out = args.out or str(trace_path.with_suffix(".png"))
        fig.savefig(out, dpi=args.dpi, bbox_inches="tight")
        print(f"Saved: {out}")
        return

    # Directory mode
    if not trace_path.is_dir():
        print(f"Not a file or directory: {trace_path}", file=sys.stderr)
        sys.exit(1)

    # Comparison mode: two groups
    if args.reads and args.reads2:
        plot_comparison(
            trace_path, args.reads, args.reads2,
            label_a=args.label or "Group A",
            label_b=args.label2 or "Group B",
            chunk=args.chunk, out=args.out, dpi=args.dpi,
        )
        return

    # Single group mode
    if args.reads:
        plot_reads(
            trace_path, args.reads,
            chunk=args.chunk, out=args.out, dpi=args.dpi,
            label=args.label,
        )
        return

    # --all mode: one PNG per read
    if args.all:
        all_files = find_trace_files(trace_path)
        print(f"Found {len(all_files)} reads in {trace_path}")
        for rid, chunks in sorted(all_files.items()):
            out_path = trace_path / f"{rid}_cloud.png"
            plot_reads(trace_path, [rid], chunk=args.chunk,
                       out=str(out_path), dpi=args.dpi)
        return

    # Default: list available reads
    all_files = find_trace_files(trace_path)
    if not all_files:
        print(f"No trace files found in {trace_path}")
        sys.exit(1)
    print(f"Found {len(all_files)} reads. Use --reads or --all to plot.\n")
    for rid, chunks in sorted(all_files.items()):
        chunk_nums = [c for c, _ in chunks]
        print(f"  {rid}  chunks={chunk_nums}")


if __name__ == "__main__":
    main()
