#!/usr/bin/env bash
#
# Systematic evaluation pipeline for panomap.
# Runs mapping, plotting, and evaluation for each read set.
#
# Usage: ./scripts/run_eval_pipeline.sh [options] [slug...]
#
# Outputs:
#   tests/data/HLA/eval_results/<slug>/
#     - anchors/          (anchor dump per read)
#     - heatmap.png       (combined heatmap)
#     - output.paf        (mapping results)
#     - eval_summary.txt  (human-readable eval)
#     - eval.json         (machine-readable eval)
#   tests/data/HLA/eval_results/SUMMARY.md  (combined report)
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PANOMAP_ROOT="$(dirname "$SCRIPT_DIR")"
PANOMAP_BIN="${PANOMAP_ROOT}/build/panomap"
PYTHON="${PYTHON:-python3}"

# Default paths
GRAPH="${PANOMAP_ROOT}/tests/data/graphs/drb1.gfa"
READS_DIR="${PANOMAP_ROOT}/tests/data/HLA/eval_reads"
RESULTS_DIR="${PANOMAP_ROOT}/tests/data/HLA/eval_results"

# Default options
FORCE=0
THREADS=5
MIN_LENGTH=0
EVENT_PIPELINE="standard"
SLUGS=""

usage() {
    echo "Usage: $0 [options] [slug...]"
    echo ""
    echo "Options:"
    echo "  --graph PATH         Path to graph file (default: tests/data/graphs/drb1.gfa)"
    echo "  -t, --threads N      Number of threads (default: 5)"
    echo "  -e, --event-pipeline Event pipeline: standard (default: standard)"
    echo "  -l, --min-length N   Filter reads shorter than N bp in eval (default: 0, no filter)"
    echo "  -f, --force          Overwrite existing results"
    echo "  -h, --help           Show this help"
    echo ""
    echo "If no slugs specified, processes all read sets in eval_reads/"
    echo ""
    echo "Examples:"
    echo "  $0                                    # Run all"
    echo "  $0 r9m_2k_ideal_dnorm                 # Run single"
    echo "  $0 -e standard                        # Use standard event pipeline"
    echo "  $0 -l 1000                            # Filter reads < 1000bp"
    echo "  $0 -f r9m_2k_ideal_dnorm              # Force rerun"
    exit 0
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --graph) GRAPH="$2"; shift 2 ;;
        -t|--threads) THREADS="$2"; shift 2 ;;
        -e|--event-pipeline) EVENT_PIPELINE="$2"; shift 2 ;;
        -l|--min-length) MIN_LENGTH="$2"; shift 2 ;;
        -f|--force) FORCE=1; shift ;;
        -h|--help) usage ;;
        -*) echo "Unknown option: $1" >&2; exit 1 ;;
        *) SLUGS="$SLUGS $1"; shift ;;
    esac
done

# Check panomap binary exists
if [[ ! -x "$PANOMAP_BIN" ]]; then
    echo "Error: panomap binary not found at $PANOMAP_BIN"
    echo "Run: cmake --build panomap/build"
    exit 1
fi

# Check graph exists
if [[ ! -f "$GRAPH" ]]; then
    echo "Error: Graph not found at $GRAPH"
    exit 1
fi

# Check reads directory exists
if [[ ! -d "$READS_DIR" ]]; then
    echo "Error: Eval reads not found at $READS_DIR"
    echo "Run: ./scripts/generate_eval_reads.sh"
    exit 1
fi

# Get list of slugs to process
if [[ -z "$SLUGS" ]]; then
    SLUGS=$(ls "$READS_DIR"/*.blow5 2>/dev/null | xargs -I{} basename {} .blow5 | sort)
fi

if [[ -z "$SLUGS" ]]; then
    echo "Error: No read sets found in $READS_DIR"
    exit 1
fi

echo "=========================================="
echo "Piru Evaluation Pipeline"
echo "=========================================="
echo "Graph:          $GRAPH"
echo "Reads:          $READS_DIR"
echo "Results:        $RESULTS_DIR"
echo "Threads:        $THREADS"
echo "Event pipeline: $EVENT_PIPELINE"
echo "Min length:     ${MIN_LENGTH}bp"
echo "Read sets:      $(echo $SLUGS | wc -w)"
echo ""

# Create results directory
mkdir -p "$RESULTS_DIR"

echo "=========================================="
echo "Processing read sets..."
echo "=========================================="

# Process each read set
process_slug() {
    local slug=$1
    local reads="${READS_DIR}/${slug}.blow5"
    local truth="${READS_DIR}/${slug}.paf"
    local outdir="${RESULTS_DIR}/${slug}"

    echo ""
    echo "--- $slug ---"

    # Check inputs exist
    if [[ ! -f "$reads" ]]; then
        echo "  SKIP: reads not found ($reads)"
        return 1
    fi
    if [[ ! -f "$truth" ]]; then
        echo "  SKIP: truth not found ($truth)"
        return 1
    fi

    # Check if already processed
    if [[ -f "$outdir/eval.json" && $FORCE -eq 0 ]]; then
        echo "  SKIP: already processed (use -f to force)"
        return 0
    fi

    # Determine model based on slug
    local model="r9.4"
    if [[ "$slug" == r10* ]]; then
        model="r10.4"
    fi

    mkdir -p "$outdir/anchors"

    # Run mapping with anchor dump
    echo "  Mapping (model=$model, event=$EVENT_PIPELINE)..."
    "$PANOMAP_BIN" map \
        --graph "$GRAPH" \
        --model "$model" \
        --event-pipeline "$EVENT_PIPELINE" \
        --linearizer path-walk \
        --clusterer dp-chain \
        --graph-type vg \
        "$reads" \
        -t "$THREADS" \
        -o "$outdir/output.paf" \
        --dump-anchors "$outdir/anchors" \
        2>&1 | grep -E "^\[" | head -10 || true

    # Generate heatmap
    echo "  Plotting heatmap..."
    "$PYTHON" "$SCRIPT_DIR/plot_anchor_heatmap.py" \
        --input "$outdir/anchors" \
        --combined \
        --paf "$outdir/output.paf" \
        --output "$outdir/heatmap.png" \
        --label "$slug" \
        2>&1 | grep -E "^(Saved|Found|Creating)" || true

    # Run evaluation
    echo "  Evaluating..."
    local eval_opts=""
    if [[ $MIN_LENGTH -gt 0 ]]; then
        eval_opts="-l $MIN_LENGTH"
    fi

    "$PANOMAP_BIN" eval \
        -t "$truth" \
        -c "$outdir/output.paf" \
        -f summary \
        $eval_opts \
        -o "$outdir/eval_summary.txt" \
        2>&1 | grep -E "^\[main\]" || true

    "$PANOMAP_BIN" eval \
        -t "$truth" \
        -c "$outdir/output.paf" \
        -f json \
        $eval_opts \
        -o "$outdir/eval.json" \
        2>/dev/null

    echo "  Done: $outdir"
}

PROCESSED=0
FAILED=0
for slug in $SLUGS; do
    if process_slug "$slug"; then
        ((PROCESSED++)) || true
    else
        ((FAILED++)) || true
    fi
done

echo ""
echo "=========================================="
echo "Generating summary report..."
echo "=========================================="

# Generate markdown summary
SUMMARY="${RESULTS_DIR}/SUMMARY.md"
{
    echo "# Piru Evaluation Results"
    echo ""
    echo "Generated: $(date '+%Y-%m-%d %H:%M:%S')"
    echo ""
    echo "## Summary Table"
    echo ""
    echo "| Read Set | Total | Mapped | Overlap Acc | Avg Length | Notes |"
    echo "|----------|-------|--------|-------------|------------|-------|"

    for slug in $SLUGS; do
        json="${RESULTS_DIR}/${slug}/eval.json"
        if [[ -f "$json" ]]; then
            # Extract values from JSON (simple grep/sed, no jq dependency)
            total=$(grep -o '"total": [0-9]*' "$json" | head -1 | grep -o '[0-9]*')
            mapped=$(grep -o '"mapped": [0-9]*' "$json" | head -1 | grep -o '[0-9]*')
            overlap_acc=$(grep -o '"overlap_accuracy": [0-9.]*' "$json" | grep -o '[0-9.]*')
            avg_len=$(grep -o '"read_length_avg_bp": [0-9]*' "$json" | grep -o '[0-9]*')

            # Calculate percentage
            if [[ -n "$overlap_acc" ]]; then
                overlap_pct=$(echo "$overlap_acc * 100" | bc -l 2>/dev/null | cut -d. -f1-2 || echo "$overlap_acc")
                overlap_pct="${overlap_pct}%"
            else
                overlap_pct="N/A"
            fi

            # Determine notes based on slug
            notes=""
            [[ "$slug" == *_ideal_* ]] && notes="ideal"
            [[ "$slug" == *_real_* ]] && notes="realistic noise"

            echo "| $slug | $total | $mapped | $overlap_pct | ${avg_len:-N/A} bp | $notes |"
        else
            echo "| $slug | - | - | - | - | not processed |"
        fi
    done

    echo ""
    echo "## Per-Read-Set Details"
    echo ""

    for slug in $SLUGS; do
        outdir="${RESULTS_DIR}/${slug}"
        if [[ -f "$outdir/eval_summary.txt" ]]; then
            echo "### $slug"
            echo ""
            echo '```'
            cat "$outdir/eval_summary.txt"
            echo '```'
            echo ""
            echo "**Heatmap:** [heatmap.png](./${slug}/heatmap.png)"
            echo ""
        fi
    done

} > "$SUMMARY"

echo "Summary written to: $SUMMARY"
echo ""
echo "=========================================="
echo "Done!"
echo "=========================================="
echo "Processed: $PROCESSED"
echo "Failed:    $FAILED"
echo ""
echo "View results:"
echo "  cat $SUMMARY"
echo "  ls $RESULTS_DIR/"
