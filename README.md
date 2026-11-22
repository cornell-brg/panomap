# Panomap

<p align="left">
  <em>Nanopore signal mapping to pangenome variation graphs</em>
</p>

<p align="left">
  <a href="LICENSE"><img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="License: MIT"></a>
</p>

## What is Panomap?

Panomap maps raw Nanopore signal reads (in SLOW5/BLOW5 format) to pangenome
variation graphs without basecalling. It works at the signal level and targets
GFA graphs built with tools like PGGB and minigraph-cactus.

Note: Panomap is an experimental proof-of-concept implementation and currently does
not scale to large genomes such as the human pangenome.

## Requirements

- **CMake** 3.16–3.31 (CMake 4.x works with one extra flag, see [Build options](#build-options))
- **C++20 compiler** (GCC 13+ recommended)
- **zlib**
- **oneTBB**: parallel index/map backend; uses a system install if present, otherwise auto-fetched during configure
- **libzstd** *(optional)*: only needed for zstd-compressed BLOW5 input

Developed and tested mainly on Red Hat Enterprise Linux 8 with GCC 13.3.1 and
CMake 3.26. If it doesn't build on your system, please open an issue or pull
request, we're happy to help.

## Quick Start

```bash
# Clone and build
git clone --recursive https://github.com/cornell-brg/panomap.git
cd panomap && mkdir build && cd build && cmake .. && make -j8

# Try the bundled example (SARS-CoV-2 pangenome + 20 simulated reads)
./panomap index -m r10.4 ../examples/covid/covid-pangenome.gfa -o covid.pirx
./panomap map --index covid.pirx ../examples/covid/reads.blow5 -o out.gaf

# --- with your own data ---

# Index a graph with a preset (r9, r9-fast, r9-small, r10, r10-fast, r10-large)
./panomap index -x r10 reference.gfa -o ref.pirx

# Map reads (GAF to stdout by default). Map uses the preset the index was built with.
./panomap map --index ref.pirx reads.blow5 -o out.gaf
```

## Build options

`cmake .. && make` works out of the box. The flags below cover special cases only.

| Flag | Default | Purpose |
|------|---------|---------|
| `-DPANOMAP_USE_ZSTD=AUTO\|ON\|OFF` | `AUTO` | zstd-compressed BLOW5 support. `AUTO` enables it when libzstd is found; `ON` fails configure if libzstd is missing; `OFF` disables. |
| `-DPANOMAP_FETCH_TBB=ON\|OFF` | `ON` | Auto-fetch oneTBB when no system install is found. `OFF` requires a system oneTBB. |
| `-DCMAKE_POLICY_VERSION_MINIMUM=3.5` | (none) | Required on **CMake 4.x**: the pinned oneTBB predates CMake 4 and otherwise errors with `Compatibility with CMake < 3.5 has been removed`. |

**zstd-compressed BLOW5.** BLOW5 records may be zlib- or zstd-compressed; zlib is
always supported. To read zstd-compressed files, install `libzstd-dev`
(Debian/Ubuntu) or `libzstd-devel` (RHEL/Fedora) before configuring; it is then
picked up automatically. Without it, mapping a zstd BLOW5 aborts with
`slow5lib has not been compiled with zstd support`.

## Presets

A preset sets the pore model, the seeds and the map settings together. `index -x <preset>`
stores the preset in the index, and `map` applies its map settings unless the command line
names another preset. Flags given after the preset override it.

| Preset | Chemistry | Use |
|--------|-----------|-----|
| `r9` | R9.4.1 | Genomes of bacterial size and up |
| `r9-fast` | R9.4.1 | About 3x faster mapping; best when reads get 4 or more chunks |
| `r9-small` | R9.4.1 | Small genomes, such as viruses |
| `r10` | R10.4.1 | Genomes of bacterial size |
| `r10-fast` | R10.4.1 | About 5x faster mapping at lower recall |
| `r10-large` | R10.4.1 | Large genomes, such as human |

## Output

`map` writes PAF for a `.paf` output name and GAF otherwise. Both formats give the query
length in signal samples (column 2) and the query start and end as event positions
(columns 3 and 4). PAF gives the chosen path as the target. In GAF, column 6 is the walk
of every path step the chain covers, columns 7 to 9 are the walk length and the chain
start and end on that walk, and `pn:Z` names the chosen path. With `map
--gaf-anchor-nodes`, column 6 lists only the nodes the anchors fall on, and columns 7 to
9 hold the path length and the positions on the path, as in PAF. Unmapped reads get a
line with `tp:A:U`.

| Tag | Type | Meaning |
|-----|------|---------|
| `tp` | A | `P` primary, `S` secondary (GAF only), `U` unmapped |
| `pn` | Z | Chosen path (GAF only; `*` if unmapped) |
| `s1` | i | Chain score |
| `an` | i | Anchors in the chain |
| `se` | f | Chain score per query event |
| `ad` | f | Anchors per reference base |
| `ws` | f | Weighted standout of the best chain |
| `nc` | i | Chains kept for the read |
| `df` | A | Decision: `G` single-chain gate, `S` standout, `F` fallback, `U` unmapped |
| `mt` | f | Map time of the read, in ms |
| `ci` | i | Chunks processed |
| `ls`, `le`, `lc` | f, f, i | 1D layout start, end and component (GAF, when the index has a layout) |

`mt` and `ci` use the names that RawHash2 and its scripts read. Read tags (`ws` to `ci`)
appear once per read, on the primary or unmapped line.

## Reproducible indexes

Two index builds of the same graph at more than one thread differ only in the 1D
layout. The layout uses path-guided SGD with lock-free parallel updates (Hogwild), as in
odgi, so its coordinates change slightly from run to run. Mapping uses the layout to
merge duplicate chains, so a few reads can map differently. `panomap inspect --compare
a.pirx b.pirx` shows which sections of two indexes differ. Builds at one thread are
identical.

## Subcommands

| Command | Description |
|---------|-------------|
| `panomap index` | Build .pirx index from GFA graph + pore model |
| `panomap map` | Map BLOW5/SLOW5 reads against index, output GAF/PAF |

Run `panomap <command> --help` for full options.

## Simulation Workflow

Generate test reads using [squigulator](https://github.com/hasindu2008/squigulator):

```bash
squigulator reference.fa -x dna-r9-min \
  -o reads.blow5 -n 20 -r 8000 \
  --sample-rate 4000 \
  -q reads.fasta -c reads.paf \
  --seed 123

./panomap index -m r9.4 reference.gfa -o ref.pirx
./panomap map --index ref.pirx reads.blow5 -o out.gaf
```

## Build & Test

Tests use [doctest](https://github.com/doctest/doctest) (fetched automatically
during configure). Build and run them with:

```bash
mkdir build && cd build
cmake ..
make -j8
make test
```

## Reproducing paper results

The [`reproduce/`](reproduce/) directory has self-contained scripts to reproduce
the panomap results in the paper from the data on Zenodo
([10.5281/zenodo.21420009](https://doi.org/10.5281/zenodo.21420009)). See
[`reproduce/README.md`](reproduce/README.md).

## Citing Panomap

If you use Panomap in your work, please cite:

```bibtex
@article{shih2026panomap,
  title   = {Panomap: Unbiased Nanopore Signal Mapping with Pangenome Variation Graphs},
  author  = {Shih, Po Jui and Sanghani, Zephan and Guarracino, Andrea and Gamaarachchi, Hasindu and Batten, Christopher},
  journal = {bioRxiv},
  year    = {2026},
  doi     = {10.64898/2026.07.10.737796},
  url     = {https://doi.org/10.64898/2026.07.10.737796}
}
```

## Acknowledgements

Panomap builds on code and ideas from several projects:

- [RawHash2](https://github.com/CMU-SAFARI/RawHash) -- signal event detection
  (originally from Scrappie), signal tokenization (adaptive quantization),
  chaining (originally from minimap2), and mapping result scoring.
- [minimap2](https://github.com/lh3/minimap2) -- seed index / hash table construction.
- [odgi](https://github.com/pangenome/odgi) -- PG-SGD 1D layout.

Bundled dependencies: [slow5lib](https://github.com/hasindu2008/slow5lib) (SLOW5/BLOW5
IO) and [kmer_models](https://github.com/nanoporetech/kmer_models) (ONT pore models).
Example and test reads are simulated with
[squigulator](https://github.com/hasindu2008/squigulator).

## License

MIT
