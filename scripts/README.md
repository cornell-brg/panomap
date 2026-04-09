# PIRU Scripts

Utility scripts for debugging and analysis.

## plot_anchor_cloud.py

Visualize anchor clouds (query_pos vs ref_coord) from trace dumps.
Shows the seed hit scatter and best chain for each read. Key diagnostic
for understanding why reads map or fail to map -- true matches show
colinear diagonals, noise shows random scatter.

### Prerequisites

```bash
pip install matplotlib numpy
```

### Enable tracing

Build with trace support and run with chain stage enabled:

```bash
cmake -S repo -B repo/build -DCMAKE_CXX_FLAGS="-DPIRU_TRACE_ENABLED"
cmake --build repo/build -j8

PIRU_TRACE_STAGES=0x40 PIRU_TRACE_DIR=/tmp/trace \
    piru map --index ref.pirx reads.blow5 --no-map-filter -o out.gaf
```

### Usage

```bash
# List available reads in a trace directory:
python scripts/plot_anchor_cloud.py /tmp/trace/

# Plot specific reads:
python scripts/plot_anchor_cloud.py /tmp/trace/ --reads S1_1 S1_5 S1_12

# Compare two groups (e.g., human noise vs microbe noise):
python scripts/plot_anchor_cloud.py /tmp/trace/ \
    --reads S1_1 S1_4 S1_6 --label "human (noise floor)" \
    --reads2 M1_1 M1_2 M1_3 --label2 "microbe (noise)" \
    --out noise_comparison.png

# Plot all reads (one PNG each):
python scripts/plot_anchor_cloud.py /tmp/trace/ --all

# Plot a specific chunk (default: last/highest):
python scripts/plot_anchor_cloud.py /tmp/trace/ --reads S1_1 --chunk 3

# Plot a single trace file directly:
python scripts/plot_anchor_cloud.py /tmp/trace/6_anchors_S1_1_chunk4
```

### Output

Each panel shows:
- Gray dots: full anchor cloud (all seed hits)
- Colored dots + line: best chain anchors
- Title: read ID, chunk number, chain score, anchor count, ws (if available)

Comparison mode (--reads + --reads2) shows two rows for side-by-side analysis.

---

## plot_signal_alignment.py

Visualization tool for debugging signal processing and mapping.

### Purpose

Helps diagnose signal processing issues by plotting:
- Node signals at different pipeline stages (raw -> tokenized -> alignment quantized)
- Read signals from BLOW5/SLOW5 files
- Anchor matches overlaid to verify signal similarity

### Prerequisites

```bash
pip install matplotlib numpy pyslow5
```

### Enable Graph Dumping

PIRU must be built with `PIRU_DUMP_GRAPHS=ON` to generate GFA annotations:

```bash
cd build
cmake -DPIRU_DUMP_GRAPHS=ON ..
cmake --build .
```

When enabled, indexing and mapping will output GFA files to `build/graph_dumps/`:
- `imported_graph.gfa` - Original graph structure
- `transformed_graph.gfa` - After VG transform (if applicable)
- `raw_signals.gfa` - Node signals after squigglization
- `tokenized.gfa` - Quantized tokens (seed extraction)
- `aln_quantized.gfa` - Alignment quantized signals

### Usage

**1. Plot node signal processing stages:**

```bash
python scripts/plot_signal_alignment.py \
    --node-id 31 \
    --graph-dir build/graph_dumps
```

Shows:
- Node sequence
- Raw signal (picoamps)
- Quantized tokens
- Alignment quantized signal

**2. Plot read signal:**

```bash
python scripts/plot_signal_alignment.py \
    --node-id 31 \
    --read-id "S1_1!gi|568815592:32578768-32589835!54!1312!-" \
    --reads tests/data/HLA/test_reads/quick_r9_1k.blow5 \
    --graph-dir build/graph_dumps
```

Shows:
- Read raw ADC signal
- Converted to picoamps

**3. Compare node vs read at anchor position:**

```bash
python scripts/plot_signal_alignment.py \
    --node-id 31 \
    --read-id "S1_1!gi|568815592:32578768-32589835!54!1312!-" \
    --reads tests/data/HLA/test_reads/quick_r9_1k.blow5 \
    --graph-dir build/graph_dumps \
    --anchor-query 482 \
    --anchor-ref 1173 \
    --window-size 100 \
    --output anchor_comparison.png
```

Shows:
- Full node and read signals
- Overlaid normalized signals at anchor window
- Correlation metric
- Quantized tokens

### Debugging Workflow

1. **Identify suspicious anchor** from mapping output:
   ```
   S1_1!...  anchors_detail=[q=482,r=1173,l=6,p=0]
   ```

2. **Find corresponding node ID** for ref_coord=1173 on path=0:
   - Check linearization coordinates (future: add query tool)
   - Or inspect GFA dumps directly

3. **Visualize anchor match**:
   ```bash
   python scripts/plot_signal_alignment.py \
       --node-id <node_id> \
       --read-id "S1_1!..." \
       --reads <reads.blow5> \
       --graph-dir build/graph_dumps \
       --anchor-query 482 \
       --anchor-ref 1173
   ```

4. **Check for issues**:
   - Low correlation (<0.5): anchor match is spurious
   - Identical tokens across all nodes: quantization collapse
   - Flat signals: event detection or normalization broken
   - Mismatched signal scales: calibration issue

### Common Issues

**"pyslow5 not available"**
```bash
pip install pyslow5
```

**"Could not load read"**
- Check read ID matches exactly (case-sensitive, special chars)
- Verify BLOW5 file path is correct

**"No GFA dumps found"**
- Rebuild with `cmake -DPIRU_DUMP_GRAPHS=ON`
- Run index or map command to generate dumps
- Check `build/graph_dumps/` directory exists

### Output Formats

- Interactive: omit `--output` to show plots with matplotlib GUI
- PNG: `--output plot.png`
- Use `.pdf` or `.svg` for vector output

---

## eval_canonical.py

Evaluate mapping accuracy using canonical 1D coordinates. Compares ground truth
(from squigulator read names) against piru's canonical interval GAF tags.

### Prerequisites

Requires piru built with component-aware dedup (dev-94). GAF output must have
`ci:f:`, `ce:f:`, `cc:i:` tags (automatic when index has 1D coords + component IDs).

### Setup

```bash
# 1. Index
piru index graph.gfa -o graph.pirx

# 2. Dump node canonical coords
piru inspect graph.pirx --dump-path-coords node_coords.tsv

# 3. Map (squigulator-simulated reads)
piru map --index graph.pirx reads.blow5 -o out.gaf
```

### Usage

```bash
python3 scripts/eval_canonical.py graph.gfa node_coords.tsv out.gaf
```

### How it works

1. Parses truth from squigulator read names: `S1_N!contig!start!end!strand`
2. Walks the GFA path for that contig in base space to find which nodes the
   truth interval covers
3. Looks up each node's canonical start/end from the coords TSV, interpolates
   to get the truth canonical interval
4. Compares against piru's `ci:f:`/`ce:f:` tags from GAF
5. Reports overlap (correct) or no overlap (wrong)

### Output

```
Total: 600
Correct: 599 (99.8%)
Wrong: 0
Unmapped: 1
No truth: 0

UNMAPPED   S1_224!gi|568815567:3779003-3792415!9413!13413!-
```

Wrong entries include truth/piru canonical intervals and node walks for
manual inspection.

### Limitations

- Only evaluates primary alignments
- Read names must follow squigulator format (`ID!contig!start!end!strand`)
- Reads with non-matching name format are silently skipped
