# PANOMAP Index Format (.pirx)

**Version**: 1.0 (pre-release)

## Overview

PANOMAP indexes are stored as a single binary `.pirx` file containing graph
topology, linearization coordinates, seed hash table, and optional 1D
canonical coordinates.

```
panomap index reference.gfa -o ref.pirx
panomap inspect ref.pirx
```

## Binary Layout

All multi-byte integers are little-endian.

```
[1. Header]
[2. Metadata]
[3. Graph nodes]
[4. Graph edges]
[5. Graph paths]
[6. Linearization coordinates]
[7. Seed store (bucket-native)]
[8. 1D canonical coordinates (optional)]
```

### 1. Header

| Offset | Size | Type | Field | Description |
|--------|------|------|-------|-------------|
| 0 | 4 | char[4] | magic | `PIRX` |
| 4 | 4 | uint32 | version | `(major << 16) \| minor` |
| 8 | 4 | uint32 | flags | Reserved (0) |

### 2. Metadata

| Type | Field | Description |
|------|-------|-------------|
| string | version | Piru version that built this index |
| uint64 | build_timestamp | Unix seconds |
| string | model_name | Pore model (e.g. `r9.4_450bps`) |
| uint32 | pore_k | Pore model k-mer size |
| string | tokenizer | Tokenizer backend (e.g. `rh2`, `landmark`) |

Strings are length-prefixed: `uint32 len` followed by `len` bytes (no null terminator).

### 3. Graph Nodes

| Type | Field | Description |
|------|-------|-------------|
| uint64 | node_count | Number of nodes |

Per node (repeated `node_count` times):

| Type | Field | Description |
|------|-------|-------------|
| string | name | Node name/ID |
| uint8 | is_reverse | 0=forward, 1=reverse |

### 4. Graph Edges

| Type | Field | Description |
|------|-------|-------------|
| uint64 | edge_count | Number of edges |

Per edge (repeated `edge_count` times):

| Type | Field | Description |
|------|-------|-------------|
| uint64 | from | Source node ID |
| uint64 | to | Target node ID |
| uint64 | overlap | Overlap in bases (0 if not tracked) |

### 5. Graph Paths

| Type | Field | Description |
|------|-------|-------------|
| uint64 | path_count | Number of paths |

Per path:

| Type | Field | Description |
|------|-------|-------------|
| string | name | Path name |
| uint64 | length | Path length in bp |
| uint64 | step_count | Number of steps |
| uint64[] | steps | Node IDs (repeated `step_count` times) |

### 6. Linearization Coordinates

| Type | Field | Description |
|------|-------|-------------|
| uint64 | node_count | Number of nodes |

Per node (repeated `node_count` times):

| Type | Field | Description |
|------|-------|-------------|
| uint64 | coord_count | Number of (path_id, ref_coord) pairs |
| (uint64, int64)[] | coords | Pairs repeated `coord_count` times |

### 7. Seed Store (Bucket-Native)

| Type | Field | Description |
|------|-------|-------------|
| string | extractor_name | Seed extractor (e.g. `kmer`) |
| uint64 | param_count | Number of key-value params |
| (string, string)[] | params | Key-value pairs (k, stride, qbits, etc.) |
| double | filter_fraction | Keep fraction (legacy, 1.0) |
| uint64 | max_hash_frequency | Max hits for any single hash |
| uint64 | freq_threshold | Frequency filter threshold |
| uint32 | bucket_bits | Log2 of number of buckets |
| uint64 | num_buckets | Total buckets (2^bucket_bits) |

Per bucket (repeated `num_buckets` times):

| Type | Field | Description |
|------|-------|-------------|
| uint32 | n_keys | Unique hashes in this bucket |
| uint32 | n_entries | Total hits in this bucket |
| uint64[] | keys | Sorted hash values (n_keys) |
| uint32[] | counts | Hit count per hash (n_keys) |
| uint32[] | offsets | Entry offset per hash (n_keys) |
| SeedEntry[] | entries | Flat hit array (n_entries) |

SeedEntry = `(uint32 node_id, uint32 offset)` = 8 bytes.

Lookup: `bucket = hash & (num_buckets - 1)`, binary search `keys[]`,
read `entries[offset..offset+count]`.

### 8. 1D Canonical Coordinates (Optional)

Present if the stream has more data after section 7.

| Type | Field | Description |
|------|-------|-------------|
| uint64 | n_coords | Number of node coordinates |
| float[] | coords | 1D positions (n_coords floats) |

## Versioning

- `major << 16 | minor` encoding
- Major bump = breaking layout change
- Minor bump = backward-compatible addition
- Current: 1.0 (pre-release, format may change without notice)

## Inspect

```
$ panomap inspect ref.pirx
index: ref.pirx
panomap_version: 0.0.1
built: 2026-04-06 20:33:29
model: r9.4_450bps
pore_k: 6
tokenizer: landmark
seed_extractor: kmer
seed_k: 4
seed_qbits: 4
...
```
