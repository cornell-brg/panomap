#!/usr/bin/env bash
#
# Generate evaluation read sets for piru testing.
# Requires: module load squigulator
#
# Usage: ./scripts/generate_eval_reads.sh [core|all]
#   core (default): 4 read sets (R9/R10 x ideal/real, 2k length)
#   all: All tiers (27 read sets)
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIRU_ROOT="$(dirname "$SCRIPT_DIR")"
REF="${PIRU_ROOT}/tests/data/HLA/DRB1-3123.fa"
OUTDIR="${PIRU_ROOT}/tests/data/HLA/eval_reads"

# Check squigulator is available
if ! command -v squigulator &> /dev/null; then
    echo "Error: squigulator not found. Run 'module load squigulator' first." >&2
    exit 1
fi

mkdir -p "$OUTDIR"

# Generate a read set
# Args: slug chemistry length noise_flags dwell_flags seed
generate() {
    local slug=$1
    local chemistry=$2
    local length=$3
    local noise_flags=$4
    local dwell_flags=$5
    local seed=$6

    local outbase="${OUTDIR}/${slug}"

    if [[ -f "${outbase}.blow5" ]]; then
        echo "  Skipping ${slug} (already exists)"
        return
    fi

    echo "  Generating ${slug}..."
    squigulator "$REF" -x "$chemistry" \
        -o "${outbase}.blow5" -n 5 -r "$length" \
        -q "${outbase}.fasta" -c "${outbase}.paf" \
        $noise_flags $dwell_flags --seed "$seed" \
        2>&1 | grep -E "^\[INFO\]|^\[main\]" || true
}

# Core: R9/R10 x ideal/real, 2k length (4 read sets)
generate_core() {
    echo "Core: R9/R10 x ideal/real, 2k length"
    generate "r9m_2k_ideal_dnorm"  "dna-r9-min"  2000 "--ideal" "" 100
    generate "r10m_2k_ideal_dnorm" "dna-r10-min" 2000 "--ideal" "" 101
    generate "r9m_2k_real_dnorm"   "dna-r9-min"  2000 ""        "" 102
    generate "r10m_2k_real_dnorm"  "dna-r10-min" 2000 ""        "" 103
}

# Extended tiers (only run with "all")
generate_extended() {
    echo "Extended: Read Length Sensitivity"
    generate "r9m_500_ideal_dnorm"  "dna-r9-min"  500   "--ideal" "" 110
    generate "r9m_1k_ideal_dnorm"   "dna-r9-min"  1000  "--ideal" "" 111
    generate "r9m_4k_ideal_dnorm"   "dna-r9-min"  4000  "--ideal" "" 112
    generate "r9m_10k_ideal_dnorm"  "dna-r9-min"  10000 "--ideal" "" 113
    generate "r9m_500_real_dnorm"   "dna-r9-min"  500   ""        "" 114
    generate "r9m_10k_real_dnorm"   "dna-r9-min"  10000 ""        "" 115
    generate "r10m_1k_ideal_dnorm"  "dna-r10-min" 1000  "--ideal" "" 116
    generate "r10m_4k_ideal_dnorm"  "dna-r10-min" 4000  "--ideal" "" 117
    generate "r10m_1k_real_dnorm"   "dna-r10-min" 1000  ""        "" 118
    generate "r10m_4k_real_dnorm"   "dna-r10-min" 4000  ""        "" 119

    echo "Extended: Noise Sensitivity"
    generate "r9m_2k_amplo_dnorm"  "dna-r9-min"  2000 "--amp-noise 0.5" "" 120
    generate "r9m_2k_amphi_dnorm"  "dna-r9-min"  2000 "--amp-noise 2.0" "" 121
    generate "r10m_2k_amplo_dnorm" "dna-r10-min" 2000 "--amp-noise 0.5" "" 122
    generate "r10m_2k_amphi_dnorm" "dna-r10-min" 2000 "--amp-noise 2.0" "" 123

    echo "Extended: Dwell/Timing Sensitivity"
    generate "r9m_2k_ideal_dfast" "dna-r9-min" 2000 "--ideal" "--dwell-mean 6.0"  130
    generate "r9m_2k_ideal_dslow" "dna-r9-min" 2000 "--ideal" "--dwell-mean 12.0" 131
    generate "r9m_2k_ideal_dvar"  "dna-r9-min" 2000 "--ideal" "--dwell-std 8.0"   132
    generate "r9m_2k_real_dfast"  "dna-r9-min" 2000 ""        "--dwell-mean 6.0"  133
    generate "r9m_2k_real_dvar"   "dna-r9-min" 2000 ""        "--dwell-std 8.0"   134

    echo "Extended: Prom vs Min Chemistry"
    generate "r9p_2k_ideal_dnorm"  "dna-r9-prom"  2000 "--ideal" "" 140
    generate "r9p_2k_real_dnorm"   "dna-r9-prom"  2000 ""        "" 141
    generate "r10p_2k_ideal_dnorm" "dna-r10-prom" 2000 "--ideal" "" 142
    generate "r10p_2k_real_dnorm"  "dna-r10-prom" 2000 ""        "" 143
}

# Main
MODE="${1:-core}"

case "$MODE" in
    core)
        generate_core
        ;;
    all)
        generate_core
        generate_extended
        ;;
    *)
        echo "Usage: $0 [core|all]" >&2
        echo "  core (default): 4 read sets (R9/R10 x ideal/real)" >&2
        echo "  all: All tiers (27 read sets)" >&2
        exit 1
        ;;
esac

echo ""
echo "Done. Generated reads in: $OUTDIR"
echo "Total files: $(ls -1 "$OUTDIR"/*.blow5 2>/dev/null | wc -l) read sets"
