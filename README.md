# PIRU

## Overview

PIRU is an experimental tool for Nanopore adaptive sampling analysis that performs direct squiggle-to-pangenome graph mapping. Rather than basecalling reads before alignment, PIRU maps raw signal data directly to graph-based reference structures, enabling real-time decision-making during sequencing.

The tool is built in modern C++ with a modular architecture that supports multiple graph formats (GFA, vg), signal file formats (slow5/blow5), and pore models (R9.4, R10.4), making it suitable for rapid prototyping and experimentation with different mapping strategies.

## Quickstart

### Dependencies

- **Build tools**: CMake ≥ 3.16 and a C++20 compiler (GCC/Clang/MSVC)
- **Core**: `zlib` development headers
- **Submodules**: `slow5lib` (required) and `libvgio` (for vg support). Run `git submodule update --init --recursive` after cloning
- **Concurrency**: oneTBB (`libtbb-dev` on Debian/Ubuntu) when `PIRU_USE_TBB=ON` (default). If missing, `PIRU_FETCH_TBB=ON` will fetch it
- **VG support**: Requires `pkg-config`, `jansson`, and `htslib` dev packages (`libjansson-dev` and `libhts-dev` on Debian/Ubuntu). Set `PIRU_ENABLE_LIBVGIO=OFF` to skip
- **Tests**: doctest is fetched automatically unless `PIRU_FETCH_DOCTEST=OFF`

### Build

```bash
# Initialize submodules (required for slow5 and libvgio)
git submodule update --init --recursive

# Configure (out-of-source) and build
cmake -S piru -B piru/build
cmake --build piru/build

# Run tests (optional)
ctest --test-dir piru/build
```

### Testing CLI Commands

To manually test the `piru index` and `piru map` commands:

1.  **Build the `piru` executable**: See the build steps above; the binary lives at `piru/build/piru`.
2.  **Create an index**: Use the `piru index` command with a test graph to generate an index directory.
    ```bash
    ./piru/build/piru index --graph=dbg --graph-k=15 \
      --output ./my_test_index piru/tests/data/graphs/sample.gfa
    ```
    This command will create a directory named `./my_test_index` containing the generated `.graph`, `.signals`, and `.seeds` files.
3.  **Run `piru map` with the new index**: Use the `piru map` command, pointing it to the index you just created and some sample reads.
    ```bash
    ./piru/build/piru map --index ./my_test_index piru/tests/data/reads/sample.blow5
    ```
    If successful, you should see the index load plus per-read seed counts printed to stdout.

**CMake options (all default to ON):**
- `PIRU_USE_TBB`: enable oneTBB backend for parallelism
- `PIRU_FETCH_TBB`: fetch oneTBB if not found on system
- `PIRU_ENABLE_LIBVGIO`: enable vg graph format support (requires jansson + htslib)
- `BUILD_TESTING`: build tests (standard CMake option)

**Common configurations:**
```bash
# Disable VG support if jansson/htslib unavailable
cmake .. -DPIRU_ENABLE_LIBVGIO=OFF

# Use system TBB instead of fetching
cmake .. -DPIRU_FETCH_TBB=OFF
```

### Basic Mapping Example

```bash
# 1. Create an index
./build/piru index --graph-k 15 --model=r10.4 --output ./my_index piru/tests/data/graphs/sample.gfa

# 2. Map reads using the created index
./build/piru map --index ./my_index piru/tests/data/reads/sample.blow5

# Show help for any command
./build/piru index --help
```

## Tool Usage

### Main Tool

```bash
./piru/build/piru [--help] [--version] <command> [options]
```

Run `piru --help` to see available subcommands and `piru <command> --help` for command-specific options.

All subcommands that support profiling (`index`, `map`, `mt-test`) accept a `-p/--profile` flag that prints timing breakdowns to stderr.

### `index` - Build Index

Load a graph file, transform it into an alignment-ready graph, squigglize/quantize signals, build seeds, and serialize the full index.

```bash
piru index [options] <graph-file>
```

**Options:**
- `-h, --help` - Show help
- `-g, --graph TYPE` - Graph type: `dbg` (default) or `vg`
- `-k, --graph-k N` - DBG k-mer size (auto-detected from overlap when possible)
- `-m, --model MODEL` - Pore model: `r9.4`, `r10.4` (builtin), or path to model file
- `-o, --output DIR` - Output directory (default `<graph-file>.piru`)
- `-t, --threads N` - Worker threads (reserved for future indexing)
- `-p, --profile` - Print timing breakdown
- `--seed-k N`, `--seed-stride N`, `--seed-filter F` - Seed extraction controls
- `--aq-backend BACKEND`, `--aq-scale SCALE` - Alignment quantization backend/override

**Examples:**
```bash
# Index a VG variation graph (e.g., HLA pangenome)
./piru/build/piru index --graph=vg --model=r9.4 \
  --output ./drb1_index piru/tests/data/graphs/drb1.vg

# Index a GFA de Bruijn graph
./piru/build/piru index --graph=dbg --graph-k=15 \
  --output ./sample_index piru/tests/data/graphs/sample.gfa

# Use custom pore model
./piru/build/piru index --graph=vg --model=/path/to/custom.model \
  --output ./my_index reference.vg
```

The command transforms the imported graph (DBG transform or VG path-guided transform), pseudo-linearizes it (chains, SCCs, superbubbles), squigglizes/quantizes signals, builds a seed store, and writes `<name>.graph`, `<name>.signals`, and `<name>.seeds` inside the output directory.

**Debugging/inspection:** Configure with `-DPIRU_DUMP_GRAPHS=ON` to emit GFA snapshots (`imported_graph.gfa`, `transformed_graph.gfa`, `raw_signals.gfa`, `aln_quantized.gfa`, `fuzzy_quantized.gfa`) for manual verification. Signal values are encoded as comma-separated sequences with tags describing the dump type and quantizer settings.

**VG Graph Indexing Details:**
- VG graphs use **path-guided transformation** with forward and reverse walks; every original edge produces forward and reverse edges.
- Nodes appearing in multiple paths with different contexts are duplicated; typical expansion is ~2x (e.g., `drb1.vg` 5111 nodes → 10222 directional nodes with k-1 overlaps).
- k-1 context is appended during path walking (both strands) so reverse nodes carry proper squiggle context.
- N bases are handled via sentinel values and are skipped during seed extraction.

### `map` - Map Reads

Load an index, stream reads, and run the read-side signal pipeline (event detection, normalization, quantization, seeding). Alignment to the graph is not wired up yet; current output reports per-read seed counts.

```bash
piru map [options] <reads-path>
```

**Options:**
- `-h, --help` - Show help
- `-i, --index DIR` - Path to index directory (required)
- `-t, --threads N` - Worker threads (reserved for future mapping)
- `-p, --profile` - Print timing breakdown

**Arguments:**
- `<reads-path>` - Input .slow5/.blow5 file or directory

**Examples:**
```bash
# Process a single blow5 file
./build/piru map reads.blow5

# Process all blow5 files in a directory
./build/piru map --profile data/reads/
```

The command scans for slow5/blow5 files and prints read metadata (ID, length, sampling rate, digitization parameters, etc.).

### `eval` - Evaluate Results

Evaluate theoretical enrichment or depletion metrics. Currently a stub awaiting implementation.

```bash
piru eval [options]
```

**Options:**
- `-h, --help` - Show help
- `-t, --truth FILE` - Ground truth/reference annotations
- `-c, --calls FILE` - Calls to evaluate

### `mt-test` - Concurrency Test

Test the parallel execution backend by spawning concurrent tasks. Useful for validating threading behavior and profiling overhead.

```bash
piru mt-test [options]
```

**Options:**
- `-h, --help` - Show help
- `-t, --threads N` - Number of parallel tasks
- `-n, --size N` - Problem size (vector length for demo)
- `-p, --profile` - Print timing breakdown

**Example:**
```bash
# Run parallel vector addition with 4 threads
./build/piru mt-test -t 4 -n 1000000 -p
```

## Developers

### Architecture

PIRU follows a modular, interface-based design that emphasizes swappable implementations and minimal dependencies:

- **Commands** (`src/commands/`) - Subcommand handlers that parse options and orchestrate pipelines
- **I/O Interfaces** (`include/io/`, `src/io/`) - Factory-based backends for reads, graphs, models, and results
  - **Reads**: `ReadProvider` with slow5/blow5 support (POD5 planned)
  - **Graphs**: `GraphLoader` for GFA and vg formats (GBZ planned)
  - **Models**: `KmerModel` with builtin R9.4/R10.4 models and file loader
  - **Results**: `ResultWriter` for GAF, GAM, JSON output (ready for integration)
- **Utilities** (`include/util/`, `src/util/`) - Concurrency (oneTBB/serial), timing, logging, metrics, signal handling
- **CLI** (`include/cli/`, `src/cli/`) - Lightweight option parsing without heavy dependencies

The codebase prefers well-maintained third-party libraries over custom implementations and maintains clear interfaces to allow experimentation with different algorithms.

### Documentation

- [ARCHITECTURE.md](docs/ARCHITECTURE.md) - High-level system architecture and module responsibilities
- [timing.md](docs/timing.md) - Profiling and timing utilities
- [threading.md](docs/threading.md) - Concurrency abstractions and executor design
- [testing.md](docs/testing.md) - Testing strategy and test organization

### Contributing

- Run `make -j $(nproc) && make check` from the build directory before submitting changes
- Follow existing code style (see `.clang-format`)
- Keep abstractions minimal and focused on current needs
- Add tests for new functionality in `tests/`
- Prefer well-maintained libraries over bespoke code

CI runs on pushes and pull requests to validate builds and tests across platforms.

## License

MIT
