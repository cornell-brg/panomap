# Reproducing panomap results from the Zenodo data

These scripts reproduce the **panomap** numbers in the paper from the Zenodo
data tarballs. Each driver is self-contained: it takes an extracted dataset
directory and the `panomap` binary, runs `index` + `map`, and evaluates.

## Prerequisites

- **Build panomap** (see the repo README) to get the `panomap` binary.
- **Python 3** -- the eval scripts use only the standard library.
- **Download + extract** the data tarballs from Zenodo
  (DOI: [10.5281/zenodo.21420009](https://doi.org/10.5281/zenodo.21420009)).
  Each `covid.tar.gz` / `yeast.tar.gz` / `zymo.tar.gz` / `hla-drb1.tar.gz`
  extracts to a `<dataset>/` folder with `graphs/ refs/ reads/`.

## Layout

```
run_covid.sh  run_yeast.sh  run_zymo.sh  run_hla-drb1.sh   # per-dataset drivers
eval/
  eval_binary.py   # covid, yeast   -> TP/FP/FN/TN, precision, recall, specificity, F1
  eval_zymo.py     # zymo           -> per-species P/R/F1 + macro_F1
  eval_drb1_recall.py  # hla-drb1   -> per-(N,chunks) recall
```

## Usage

Extract a dataset, then run its driver pointed at the extracted directory:

```bash
tar xzf covid.tar.gz
./run_covid.sh    covid     /path/to/panomap
./run_yeast.sh    yeast     /path/to/panomap
./run_zymo.sh     zymo      /path/to/panomap
./run_hla-drb1.sh hla-drb1  /path/to/panomap 07     # tier: 15 (near) | 12 (moderate) | 07 (far)
```

Each writes to `./build/<dataset>/` and prints a metrics table.

Optional env overrides:
- `BUILD=pggb|mc`  -- graph builder (covid/yeast/zymo; default `pggb`)
- `N_VALUES="1 2 4 8"`  -- pangenome sizes (default all)
- `CHUNKS="1 2 4 8"`  -- max-chunks grid (default all)
- `PY=python3`  -- python for the eval scripts

## Metrics reported

- **covid, yeast** (binary): TP/FP/FN/TN, precision, recall, specificity, F1 per (N, chunks).
- **zymo** (multi-class): per-species P/R/F1 and macro_F1 per (N, chunks), plus confusion matrices.
- **hla-drb1** (recall): per-(N, chunks) recall, for the chosen divergence tier.

Per-dataset map arguments are hardcoded in each driver (from the paper
Supplementary S4).
