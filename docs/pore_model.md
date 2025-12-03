# Pore Model

PIRU uses a simple lookup-table pore model that maps k-mers to expected nanopore current means (“squiggle” values). The model API is intentionally small: given a k-mer string, return its expected mean current. Standard deviation is not tracked at this stage.

## Built-in models
- `r9.4` / `r9.4_450bps`: 6-mer DNA model.
- `r10.4` / `r10.4_400bps`: 9-mer DNA model.

These are checked in as header-only tables under `src/io/models/builtin_model_*.hpp` and are always available without extra files.

## Loading external models
- `load_model_from_file(path)` accepts text tables with either:
  - Legacy R9 format: header row with `kmer level_mean ...` followed by rows of k-mer and mean values; extra columns are ignored.
  - Two-column format: `kmer <mean>` per line.
- K-mer length is inferred from the first parsed k-mer; the whole table must match that k.
- Unknown or unreadable files log an error and return `nullptr`.

### File examples
Legacy R9 (header + many columns):
```
kmer    level_mean  level_stdv  sd_mean  sd_stdv  weight
AAAAAA  86.486336   1.517846   0.94     0.61     4739.55
AAAACC  75.706298   1.705015   1.00     0.70     3403.70
```

Two-column (R10-style):
```
AAAAAAAAA   -1.8424464
AAAAAAAAC   -1.6519798
```

## Command-line usage
- `piru index -m <model>` accepts a built-in name (e.g., `r10.4`) or a model file path. Defaults to `r10.4`.

## Interface
`include/io/models/model.hpp`:
- `std::string name() const`
- `int k() const`
- `bool lookup(const std::string& kmer, double& mean) const`

Factory helpers in `include/io/models/model_factory.hpp`:
- `load_builtin_model(name)` for `r9.4` / `r10.4` tables.
- `load_model_from_file(path)` for runtime tables (formats above).

### Typical usage
```c++
auto model = piru::io::load_builtin_model("r10.4");
if (!model) { /* handle error */ }

const std::string kmer = "AAAAAAAAA";  // length must match model->k()
double mean = 0.0;
if (model->lookup(kmer, mean)) {
    // use mean in scaling/likelihood computations
}
```

Load from file instead of built-in:
```c++
auto model = piru::io::load_model_from_file("models/custom.model");
if (!model) { /* handle error */ }
double mean = 0.0;
if (model->lookup("AAAAAA", mean)) {
    // consume mean
}
```

## Scope notes
- Only DNA models are supported for now; RNA tables are not loaded.
- The interface returns mean-only values; downstream logic can add scaling or variability as needed.
