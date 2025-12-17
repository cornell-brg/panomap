# Results I/O Overview

This doc summarizes the alignment result formats, structures, and writers.

## Output Formats

### PAF (Pairwise mApping Format)

PAF is a TAB-delimited format for sequence-to-sequence alignments. Used when mapping to linear paths.

**12 Mandatory Columns:**

| Col | Field | Type | Description |
|-----|-------|------|-------------|
| 1 | query_name | string | Query/read name |
| 2 | query_len | int | Query sequence length |
| 3 | query_start | int | Query start (0-based, closed) |
| 4 | query_end | int | Query end (0-based, open) |
| 5 | strand | char | `+` if same strand, `-` if opposite |
| 6 | target_path | string | Target path/sequence name |
| 7 | target_len | int | Target path length |
| 8 | target_start | int | Target start (0-based) |
| 9 | target_end | int | Target end (0-based) |
| 10 | matches | int | Number of matching bases |
| 11 | block_len | int | Alignment block length |
| 12 | mapq | int | Mapping quality (0-255) |

**Example:**
```
read1	1000	50	800	+	HLA-A*01:01:01:01	3000	100	850	700	750	60	tp:A:P	cs:i:500	an:i:25
```

### GAF (Graph Alignment Format)

GAF extends PAF for graph alignments. Column 6 contains the graph traversal path instead of a sequence name.

**12 Mandatory Columns:**

| Col | Field | Type | Description |
|-----|-------|------|-------------|
| 1 | query_name | string | Query/read name |
| 2 | query_len | int | Query sequence length |
| 3 | query_start | int | Query start (0-based, closed) |
| 4 | query_end | int | Query end (0-based, open) |
| 5 | strand | char | `+` if same strand, `-` if opposite |
| 6 | graph_path | string | **Graph traversal** (e.g., `>n1>n2>n3`) |
| 7 | path_len | int | Path length through the graph |
| 8 | path_start | int | Start position on path (0-based) |
| 9 | path_end | int | End position on path (0-based) |
| 10 | matches | int | Number of matching bases |
| 11 | block_len | int | Alignment block length |
| 12 | mapq | int | Mapping quality (0-255) |

**Graph Path Format (Column 6):**
- Uses `>` for forward orientation, `<` for reverse
- Node IDs are the original graph node identifiers
- Example: `>s1>s2>s3` means forward through nodes s1, s2, s3

**Example:**
```
read1	1000	50	800	+	>s1>s2>s3	3000	100	850	700	750	60	pn:Z:HLA-A*01:01:01:01	tp:A:P	cs:i:500	an:i:25
```

### Key Difference: PAF vs GAF

| | PAF | GAF |
|--|-----|-----|
| Column 6 | Path/sequence name | Graph traversal (`>n1>n2>n3`) |
| Use case | Linear reference | Graph reference |
| Path name | In column 6 | In `pn:Z:` tag |

## Optional Tags

Tags follow SAM-style format: `TAG:TYPE:VALUE`

### Tags produced by piru

| Tag | Type | Description |
|-----|------|-------------|
| `tp` | A (char) | Alignment type: `P` = primary, `S` = secondary |
| `cs` | i (int) | Chain score from DP chaining |
| `an` | i (int) | Number of anchors in the chain |
| `pn` | Z (string) | Path name (GAF only, since col 6 is graph path) |

### Standard tags (future)

| Tag | Type | Description |
|-----|------|-------------|
| `cg` | Z (string) | CIGAR string |
| `NM` | i (int) | Edit distance |
| `AS` | i (int) | Alignment score |

## Data Model

`AlignmentResult` (`include/io/results/result.hpp`):

```cpp
struct AlignmentResult {
    // Query info
    std::string query_name;
    std::uint64_t query_length, query_start, query_end;
    std::string query_sequence, query_quality;

    char strand;  // '+' or '-'

    // Target info
    std::string target_path;   // Path name (PAF col 6)
    std::string graph_path;    // Graph traversal (GAF col 6)
    std::uint64_t target_length, target_start, target_end;

    // Alignment stats
    std::uint64_t matches, alignment_block_length;
    int mapq;

    // Detailed mappings (for GAM/JSON)
    std::vector<Mapping> mappings;

    // Optional tags as formatted strings
    std::vector<std::string> optional_fields;
};
```

## Writers

| Extension | Writer | Description |
|-----------|--------|-------------|
| `.paf` | PafWriter | TAB-delimited PAF (no dependencies) |
| `.gaf` | GafWriter | TAB-delimited GAF (no dependencies) |
| `.gam` | GamWriter | vg protobuf (requires libvgio) |
| `.json` | JsonWriter | vg JSON (requires libvgio) |

### Factory Usage

```cpp
#include "io/results/result_writer_factory.hpp"

// Auto-detect format from extension
auto writer = piru::io::make_result_writer("output.paf");

// Or with explicit format override
auto writer = piru::io::make_result_writer("output", "gaf");

piru::io::AlignmentResult r;
// ... populate r ...
writer->write(r);
```

### CLI Usage

```bash
# PAF output (path name in column 6)
./piru map --graph graph.gfa --model r9.4 -o results.paf reads.blow5

# GAF output (graph traversal in column 6)
./piru map --graph graph.gfa --model r9.4 -o results.gaf reads.blow5

# Format override
./piru map --graph graph.gfa --model r9.4 -o results --output-format paf reads.blow5
```

## References

- [PAF format specification](https://github.com/lh3/miniasm/blob/master/PAF.md)
- [GAF format specification](https://github.com/lh3/gfatools/blob/master/doc/rGFA.md)
- [minimap2 documentation](https://lh3.github.io/minimap2/minimap2.html)
