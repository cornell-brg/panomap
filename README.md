# PIRU

PIRU is an experimental squiggle-to-pangenome graph mapper for Nanopore adaptive sampling, written in modern C++ with pluggable backends.

## Status
- CLI scaffold with subcommands (`index-vg`, `index-dbg`, `map`, `eval`, `mt-test`).
- Concurrency abstraction with oneTBB backend (required when enabled) and serial fallback.
- Optional timing/profiling via start/stop API and profiling macros.
- Doctest/CTest harness with sample tests; CI builds and runs tests.

## Layout
- `src/` — commands, concurrency backend, timing, signal handlers, CLI parsing.
- `include/` — public headers.
- `tests/` — doctest-based unit tests.
- `docs/` — developer docs (timing, threading, testing).
- `plans/` — planning notes.

## Quickstart
Dependencies: CMake ≥ 3.16, C++17 compiler, and oneTBB (`libtbb-dev` on Debian/Ubuntu). Doctest is fetched automatically if not installed.

```bash
cmake -S . -B build -DBUILD_TESTING=ON
cmake --build build
ctest --test-dir build --output-on-failure   # or: cmake --build build --target check
```

Key options:
- `PIRU_USE_TBB` (default ON): enable oneTBB backend.
- `PIRU_FETCH_TBB` (default ON): fetch oneTBB if not installed.
- `PIRU_FETCH_DOCTEST` (default ON): fetch doctest if not installed.

## Profiling
- Subcommands that support profiling (`index-*`, `map`, `mt-test`) accept `-p/--profile` and print a timing summary to stderr.
- Use `PIRU_PROFILE_START/STOP` macros in code to add scoped timings; see `docs/timing.md`.

## Usage
- Top-level help:
  ```bash
  ./build/piru --help
  ```
- Version:
  ```bash
  ./build/piru --version
  ```
- Subcommand help:
  ```bash
  ./build/piru <command> --help
  ```

### Subcommands
- `index-vg` — stub for vg-based indexing. Options: `-r/--reads`, `-g/--graph`, `-t/--threads`, `-p/--profile`.
- `index-dbg` — stub for de Bruijn graph indexing. Options: `-r/--reads`, `-k/--kmer`, `-t/--threads`, `-p/--profile`.
- `map` — stub for mapping. Options: `-r/--reads`, `-i/--index`, `-t/--threads`, `-p/--profile`.
- `eval` — stub for evaluating calls. Options: `-t/--truth`, `-c/--calls`.
- `mt-test` — parallel vector addition demo. Options: `-t/--threads` (thread hint/tasks), `-n/--size` (vector length), `-p/--profile` (print timing).

### Examples
```bash
# Show top-level help
./build/piru --help

# Inspect options for a subcommand
./build/piru index-vg --help

# Run the concurrency demo with profiling
./build/piru mt-test -t 4 -n 1000000 -p
```

## Contributing
- Keep modules swappable via clear interfaces; prefer well-maintained libraries over bespoke code.
- Run `cmake --build build` and `ctest --test-dir build` before sending changes. A basic CI workflow builds and tests on pushes/PRs.
