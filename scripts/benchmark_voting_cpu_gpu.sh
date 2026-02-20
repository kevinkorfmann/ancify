#!/usr/bin/env bash
# ──────────────────────────────────────────────────────────────
# Benchmark voting run-time: CPU vs GPU for several human chromosomes.
# Runs voting only, with timing stratified by phase:
#   Phase 1 = projection (AXT → projected FASTA)
#   Phase 2 = ancestral calling (voting)
# Records phase1_sec, phase2_sec, and total time_sec per (chromosome, backend),
# then plots grouped bar charts.
#
# Env:
#   CHROMOSOMES  space-separated list (default: chr20 chr21 chr22)
#   WORK_DIR     working directory (default: repo/ancify_voting_bench)
#   ANCIFY_CPUS  workers for CPU runs (default: 4)
#   SKIP_PLOT    set to 1 to skip plotting (e.g. if no matplotlib)
#
# Usage:
#   ./scripts/benchmark_voting_cpu_gpu.sh
#   CHROMOSOMES="chr21 chr22" ./scripts/benchmark_voting_cpu_gpu.sh
# ──────────────────────────────────────────────────────────────
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

CHROMOSOMES="${CHROMOSOMES:-chr20 chr21 chr22}"
export WORK_DIR="${WORK_DIR:-$REPO_DIR/ancify_voting_bench}"
export METHOD=voting
export ANCIFY_CPUS="${ANCIFY_CPUS:-4}"
SKIP_PLOT="${SKIP_PLOT:-0}"

TIMINGS_CSV="$WORK_DIR/voting_timings.csv"
PLOT_PATH="$WORK_DIR/voting_cpu_vs_gpu.png"

echo "╔═══════════════════════════════════════════════════════════╗"
echo "║  Voting benchmark: CPU vs GPU (human hg38)                ║"
echo "║  Timing stratified: Phase 1 (project) + Phase 2 (voting)  ║"
echo "║  Chromosomes: $CHROMOSOMES"
echo "║  Work dir:    $WORK_DIR"
echo "╚═══════════════════════════════════════════════════════════╝"
echo ""

mkdir -p "$WORK_DIR"
echo "chromosome,backend,phase1_sec,phase2_sec,time_sec" > "$TIMINGS_CSV"

for CHROM in $CHROMOSOMES; do
  export CHROM
  for BACKEND in cpu gpu; do
    export BACKEND
    CONFIG="$WORK_DIR/config_${CHROM}_voting.yaml"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "  $CHROM  —  backend: $BACKEND"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    CONFIG_ONLY=1 "$SCRIPT_DIR/run_hg38.sh" || { echo "  FAILED (config)"; echo "$CHROM,$BACKEND,-1,-1,-1" >> "$TIMINGS_CSV"; continue; }
    if [[ ! -f "$CONFIG" ]]; then
      echo "  FAILED (no config)"
      echo "$CHROM,$BACKEND,-1,-1,-1" >> "$TIMINGS_CSV"
      continue
    fi
    cd "$REPO_DIR"
    # Phase 1: project
    start1=$(date +%s)
    ancify project -c "$CONFIG" || { echo "  FAILED (project)"; echo "$CHROM,$BACKEND,-1,-1,-1" >> "$TIMINGS_CSV"; continue; }
    end1=$(date +%s)
    phase1=$((end1 - start1))
    echo "  Phase 1 (project): ${phase1}s"
    # Phase 2: call (voting)
    start2=$(date +%s)
    ancify call -c "$CONFIG" || { echo "  FAILED (call)"; echo "$CHROM,$BACKEND,$phase1,-1,-1" >> "$TIMINGS_CSV"; continue; }
    end2=$(date +%s)
    phase2=$((end2 - start2))
    total=$((phase1 + phase2))
    echo "$CHROM,$BACKEND,$phase1,$phase2,$total" >> "$TIMINGS_CSV"
    echo "  Phase 2 (voting): ${phase2}s  → total ${total}s"
    echo ""
  done
done

echo "╔═══════════════════════════════════════════════════════════╗"
echo "║  Timings written to $TIMINGS_CSV"
echo "╚═══════════════════════════════════════════════════════════╝"
echo ""
cat "$TIMINGS_CSV"

if [[ "$SKIP_PLOT" == "1" ]]; then
  echo "SKIP_PLOT=1 — not generating plot."
  exit 0
fi

echo "Generating plot..."
python3 "$SCRIPT_DIR/plot_voting_benchmark.py" "$TIMINGS_CSV" "$PLOT_PATH" || {
  echo "Plot failed (install matplotlib: pip install 'ancify[evaluate]')."
  exit 0
}
echo "Plot saved: $PLOT_PATH"
