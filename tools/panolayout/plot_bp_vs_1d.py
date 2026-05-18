#!/usr/bin/env python3
"""plot_bp_vs_1d.py -- scatter (cum_bp, 1D coord) per path from panolayout's
<prefix>.path_steps.tsv. Each step contributes a single dot at the node's
canon_start (or canon_end if orient='-') so the slope reveals layout vs
traversal direction.
"""
import argparse
import csv
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from collections import defaultdict


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("path_steps_tsv")
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--species", default="all",
                    help="Filter paths by prefix (split on '#'); 'all' = no filter")
    ap.add_argument("--title", default=None)
    ap.add_argument("--alpha", type=float, default=0.4)
    ap.add_argument("--markersize", type=float, default=2.0)
    ap.add_argument("--subsample", type=int, default=1, help="Plot every Nth step (default 1)")
    ap.add_argument("--cmap", default="gist_rainbow",
                    help="Matplotlib colormap for per-path colors (default gist_rainbow)")
    ap.add_argument("--no-legend", action="store_true",
                    help="Suppress legend (useful for very many paths)")
    args = ap.parse_args()

    paths = defaultdict(lambda: ([], []))  # name -> (xs, ys)
    step_ctr = defaultdict(int)
    with open(args.path_steps_tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            name = row["path_name"]
            if args.species != "all":
                sp = name.split("#", 1)[0]
                if sp != args.species:
                    continue
            idx = step_ctr[name]
            step_ctr[name] += 1
            if args.subsample > 1 and idx % args.subsample != 0:
                continue
            x = float(row["cum_bp"])
            y = float(row["coord_start"]) if row["orient"] == "+" else float(row["coord_end"])
            paths[name][0].append(x)
            paths[name][1].append(y)

    if not paths:
        print(f"No paths matched filter '{args.species}'")
        return

    fig, ax = plt.subplots(figsize=(12, 9))
    ordered_names = sorted(paths.keys())
    n_paths = len(ordered_names)
    cmap = cm.get_cmap(args.cmap, max(n_paths, 1))
    for i, name in enumerate(ordered_names):
        xs, ys = paths[name]
        ax.scatter(xs, ys, s=args.markersize, alpha=args.alpha,
                   color=cmap(i), label=name)
    ax.set_xlabel("cumulative bp along path")
    ax.set_ylabel("PG-SGD 1D coord")
    ax.set_title(args.title or f"bp vs 1D ({args.path_steps_tsv})")
    if not args.no_legend:
        # Scale legend font + columns by path count so it fits.
        if n_paths <= 10:
            fontsize, ncol = 9, 1
        elif n_paths <= 40:
            fontsize, ncol = 6, 2
        elif n_paths <= 100:
            fontsize, ncol = 5, 3
        else:
            fontsize, ncol = 4, 4
        ax.legend(markerscale=4, fontsize=fontsize, loc="center left",
                  bbox_to_anchor=(1.01, 0.5), ncol=ncol, framealpha=0.9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(args.out, dpi=120, bbox_inches="tight")
    print(f"-> {args.out} ({len(paths)} paths)")


if __name__ == "__main__":
    main()
