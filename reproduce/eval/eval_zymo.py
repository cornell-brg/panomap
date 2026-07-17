#!/usr/bin/env python3
"""eval_zymo_multispecies.py -- multi-species classification metrics for Zymo.

Walks build_dir/results/runs/zymo-N{1,2,4,8}/chunks_{1,2,4,8}/on_target.gaf,
parses each read's `pn:Z:` (predicted path-species) and joins with the ground-
truth species TSV (sigmoni-zymo.read-species.tsv).

Emits in build_dir/results/summary/:
  per_species_metrics.tsv  -- per (N, chunks, species): TP, FP, FN, P, R, F1,
                              n_true, n_predicted, n_correct
  macro_summary.tsv        -- per (N, chunks): macro_F1, mapped_frac, n_total
  confusion_N{N}_C{C}.tsv  -- 8x8 confusion matrix per (N, chunks_max)

Path-name parsing: piru `pn:Z:` follows the format `<species>#<accession>` from
the multi-component zymo-7sp-pggb concat. Unmapped reads (`pn:Z:*`) bucket as
'unmapped'.

Usage:
  eval_zymo_multispecies.py <build_dir> --truth /path/to/species.tsv
"""

import argparse
import sys
from collections import defaultdict
from pathlib import Path

SPECIES_LIST = [
    "ecoli", "bsubtilis", "efaecalis", "lmonocytogenes",
    "paeruginosa", "senterica", "saureus",
]

# Map TSV contig name -> our species code
CONTIG_RULES = [
    ("Escherichia_coli", "ecoli"),
    ("Bacillus", "bsubtilis"),
    ("BS.", "bsubtilis"),
    ("Enterococcus", "efaecalis"),
    ("Listeria", "lmonocytogenes"),
    ("Pseudomonas", "paeruginosa"),
    ("Salmonella", "senterica"),
    ("Staphylococcus", "saureus"),
]


def contig_to_species(contig):
    for substr, sp in CONTIG_RULES:
        if substr in contig:
            return sp
    return "other"


def path_to_species(pn):
    if not pn or pn == "*":
        return None
    return pn.split("#", 1)[0]


def load_truth(path):
    """read_id -> species (TSV: read_id\tcontig\tmapq)."""
    truth = {}
    with open(path) as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 3:
                continue
            sp = contig_to_species(cols[1])
            if sp == "other":
                continue
            truth[cols[0]] = sp
    return truth


def parse_gaf_predictions(gaf_path):
    """Yield (read_id, mapped, predicted_species)."""
    with open(gaf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 12:
                continue
            qname = cols[0]
            tags = cols[12:]
            mapped = "tp:A:P" in tags
            pn = ""
            for t in tags:
                if t.startswith("pn:Z:"):
                    pn = t[5:]
                    break
            predicted = path_to_species(pn) if mapped else "unmapped"
            yield qname, mapped, predicted


def compute_metrics(gaf_path, truth):
    """Return (confusion, n_total_in_truth, n_mapped, n_unmapped)."""
    confusion = defaultdict(int)  # (true_sp, predicted_sp) -> count
    n_total = 0
    n_mapped = 0
    n_unmapped = 0
    seen_ids = set()
    for qname, mapped, predicted in parse_gaf_predictions(gaf_path):
        true_sp = truth.get(qname)
        if true_sp is None:
            continue  # not in our ground truth
        seen_ids.add(qname)
        n_total += 1
        if mapped:
            n_mapped += 1
        else:
            n_unmapped += 1
        confusion[(true_sp, predicted)] += 1
    # Reads in truth but missing from GAF: count as unmapped
    missing = set(truth) - seen_ids
    for rid in missing:
        true_sp = truth[rid]
        n_total += 1
        n_unmapped += 1
        confusion[(true_sp, "unmapped")] += 1
    return confusion, n_total, n_mapped, n_unmapped


def per_species_from_confusion(confusion):
    """Return list of dicts per species: TP, FP, FN, P, R, F1, n_true, n_pred."""
    rows = []
    for sp in SPECIES_LIST:
        TP = confusion.get((sp, sp), 0)
        FN = sum(c for (t, p), c in confusion.items() if t == sp and p != sp)
        FP = sum(c for (t, p), c in confusion.items() if t != sp and p == sp)
        n_true = TP + FN
        n_pred = TP + FP
        P = TP / n_pred if n_pred > 0 else 0.0
        R = TP / n_true if n_true > 0 else 0.0
        F1 = 2 * P * R / (P + R) if (P + R) > 0 else 0.0
        rows.append({
            "species": sp,
            "TP": TP, "FP": FP, "FN": FN,
            "n_true": n_true, "n_predicted": n_pred,
            "P": P, "R": R, "F1": F1,
        })
    return rows


def write_tsv(path, rows, fields=None):
    if not rows:
        path.write_text("")
        return
    if fields is None:
        fields = list(rows[0].keys())
    with path.open("w") as f:
        f.write("\t".join(fields) + "\n")
        for r in rows:
            vals = []
            for k in fields:
                v = r[k]
                if isinstance(v, float):
                    vals.append(f"{v:.4f}" if abs(v) < 1 else f"{v:.3f}")
                else:
                    vals.append(str(v))
            f.write("\t".join(vals) + "\n")


def write_confusion(path, confusion):
    """8x8 confusion matrix (rows = true species + 'unmapped' col)."""
    cols = SPECIES_LIST + ["unmapped"]
    with path.open("w") as f:
        f.write("true\\predicted\t" + "\t".join(cols) + "\n")
        for true_sp in SPECIES_LIST:
            row = [true_sp]
            for col_sp in cols:
                row.append(str(confusion.get((true_sp, col_sp), 0)))
            f.write("\t".join(row) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("build_dir")
    ap.add_argument("--truth", required=True,
                    help="Species TSV (read_id\\tcontig\\tmapq)")
    args = ap.parse_args()

    bd = Path(args.build_dir).resolve()
    runs_root = bd / "results" / "runs"
    summary = bd / "results" / "summary"
    summary.mkdir(parents=True, exist_ok=True)

    truth = load_truth(args.truth)
    print(f"Loaded {len(truth):,} ground-truth labels from {args.truth}")

    per_sp_rows = []
    macro_rows = []
    for n_dir in sorted(p for p in runs_root.iterdir()
                        if p.is_dir() and p.name.startswith("zymo-N")):
        try:
            n_val = int(n_dir.name.removeprefix("zymo-N"))
        except ValueError:
            continue
        chunk_dirs = sorted(
            n_dir.glob("chunks_*"),
            key=lambda p: int(p.name.removeprefix("chunks_")))
        for ck_dir in chunk_dirs:
            chunks_max = int(ck_dir.name.removeprefix("chunks_"))
            gaf = ck_dir / "on_target.gaf"
            if not gaf.exists():
                continue
            confusion, n_total, n_mapped, n_unmapped = compute_metrics(gaf, truth)
            sp_rows = per_species_from_confusion(confusion)

            # Per-species rows tagged with N + chunks
            for r in sp_rows:
                r2 = {"N": n_val, "chunks_max": chunks_max, **r}
                per_sp_rows.append(r2)

            # Macro
            f1_vals = [r["F1"] for r in sp_rows]
            macro_F1 = sum(f1_vals) / len(f1_vals) if f1_vals else 0.0
            macro_rows.append({
                "N": n_val, "chunks_max": chunks_max,
                "n_total": n_total, "n_mapped": n_mapped, "n_unmapped": n_unmapped,
                "mapped_frac": n_mapped / n_total if n_total > 0 else 0.0,
                "macro_F1": macro_F1,
            })

            # Confusion matrix per (N, chunks)
            cm_path = summary / f"confusion_N{n_val}_C{chunks_max}.tsv"
            write_confusion(cm_path, confusion)

    write_tsv(summary / "per_species_metrics.tsv", per_sp_rows,
              fields=["N", "chunks_max", "species", "TP", "FP", "FN",
                      "n_true", "n_predicted", "P", "R", "F1"])
    write_tsv(summary / "macro_summary.tsv", macro_rows,
              fields=["N", "chunks_max", "n_total", "n_mapped", "n_unmapped",
                      "mapped_frac", "macro_F1"])

    print(f"Wrote to {summary}/")
    print(f"  per_species_metrics.tsv ({len(per_sp_rows)} rows)")
    print(f"  macro_summary.tsv       ({len(macro_rows)} rows)")
    print(f"  confusion_N*_C*.tsv     (per-(N,chunks) confusion matrices)")


if __name__ == "__main__":
    main()
