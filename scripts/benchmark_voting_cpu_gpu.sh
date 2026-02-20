#!/usr/bin/env bash
# ──────────────────────────────────────────────────────────────
# Benchmark voting run-time: CPU vs GPU across human chromosomes.
#
# Uses the full 4-outgroup configuration (bonobo, chimp, gorilla,
# macaque) so the GPU has enough work to demonstrate acceleration.
# With only 2 outgroups the per-position compute is too trivial
# for GPU kernel launch overhead to pay off.
#
# Timing is stratified by phase:
#   Phase 1 = projection (AXT → projected FASTA)
#   Phase 2 = ancestral calling (voting)
#
# Env:
#   CHROMOSOMES  space-separated list
#                (default: chr1 chr2 … chr22 chrX)
#   WORK_DIR     working directory (default: repo/ancify_voting_bench)
#   ANCIFY_CPUS  workers for CPU runs (default: 4)
#   SKIP_PLOT    set to 1 to skip plotting (e.g. if no matplotlib)
#
# Usage:
#   ./scripts/benchmark_voting_cpu_gpu.sh          # all chromosomes
#   CHROMOSOMES="chr1" ./scripts/benchmark_voting_cpu_gpu.sh
# ──────────────────────────────────────────────────────────────
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

if [[ -z "${CHROMOSOMES:-}" ]]; then
  CHROMOSOMES="$(printf 'chr%s ' {1..22})chrX"
fi
WORK_DIR="${WORK_DIR:-$REPO_DIR/ancify_voting_bench}"
ANCIFY_CPUS="${ANCIFY_CPUS:-4}"
SKIP_PLOT="${SKIP_PLOT:-0}"

BASE_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38"

N_CHROMS=$(echo $CHROMOSOMES | wc -w)
TIMINGS_CSV="$WORK_DIR/voting_timings.csv"
PLOT_PATH="$WORK_DIR/voting_cpu_vs_gpu.png"

echo "╔══════════════════════════════════════════════════════════╗"
echo "║  Voting benchmark: CPU vs GPU (human hg38, 4 outgroups) ║"
echo "║  Phases: 1 (projection) + 2 (voting)                   ║"
echo "║  Chromosomes: $N_CHROMS  ($(echo $CHROMOSOMES | cut -d' ' -f1) … $(echo $CHROMOSOMES | awk '{print $NF}'))"
echo "║  CPUs:        $ANCIFY_CPUS"
echo "║  Work dir:    $WORK_DIR"
echo "╚══════════════════════════════════════════════════════════╝"
echo ""

mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

# ── Download data (all 4 outgroups) ──────────────────────────
download() {
  local file="$1" subdir="$2"
  if [[ ! -f "$file" ]]; then
    echo "[download] $file..."
    wget -q --show-progress "$BASE_URL/$subdir/$file" -O "$file" || {
      echo "  wget failed; try: curl -L $BASE_URL/$subdir/$file -o $file"
      exit 1
    }
  fi
}

if [[ ! -f hg38.chrom.sizes ]]; then
  echo "[download] hg38 chromosome sizes..."
  curl -sL "$BASE_URL/bigZips/hg38.chrom.sizes" -o hg38.chrom.sizes
fi

download hg38.panPan3.net.axt.gz   vsPanPan3
download hg38.panTro6.net.axt.gz   vsPanTro6
download hg38.gorGor6.net.axt.gz   vsGorGor6
download hg38.rheMac10.net.axt.gz  vsRheMac10

echo ""

# ── Run benchmark ────────────────────────────────────────────
echo "chromosome,backend,phase1_sec,phase2_sec,time_sec" > "$TIMINGS_CSV"

DONE=0
FAILED=0
for CHROM in $CHROMOSOMES; do
  # Write per-chromosome config with all 4 outgroups
  awk -v c="$CHROM" '$1 == c' hg38.chrom.sizes > "$WORK_DIR/chromoLens_${CHROM}.txt"
  if [[ ! -s "$WORK_DIR/chromoLens_${CHROM}.txt" ]]; then
    echo "WARNING: $CHROM not found in hg38.chrom.sizes — skipping"
    continue
  fi

  for BACKEND in cpu gpu; do
    OUT_DIR="$WORK_DIR/output_${CHROM}/voting"
    CONFIG="$WORK_DIR/config_${CHROM}_voting.yaml"

    cat > "$CONFIG" <<EOF
focal_species: human
chromosome_lengths: $WORK_DIR/chromoLens_${CHROM}.txt
chromosomes:
  - $CHROM
outgroups:
  inner:
    - name: bonobo
      alignment: $WORK_DIR/hg38.panPan3.net.axt.gz
    - name: chimp
      alignment: $WORK_DIR/hg38.panTro6.net.axt.gz
    - name: gorilla
      alignment: $WORK_DIR/hg38.gorGor6.net.axt.gz
  outer:
    - name: macaque
      alignment: $WORK_DIR/hg38.rheMac10.net.axt.gz
work_dir: $WORK_DIR
output_dir: $OUT_DIR
num_cpus: $ANCIFY_CPUS
method: voting
backend: $BACKEND
EOF

    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "  [$((DONE/2 + 1))/$N_CHROMS]  $CHROM  —  $BACKEND  (4 outgroups)"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

    cd "$REPO_DIR"
    start1=$(date +%s)
    ancify project -c "$CONFIG" || { echo "  FAILED (project)"; echo "$CHROM,$BACKEND,-1,-1,-1" >> "$TIMINGS_CSV"; FAILED=$((FAILED+1)); continue; }
    end1=$(date +%s)
    phase1=$((end1 - start1))

    start2=$(date +%s)
    ancify call -c "$CONFIG" || { echo "  FAILED (call)"; echo "$CHROM,$BACKEND,$phase1,-1,-1" >> "$TIMINGS_CSV"; FAILED=$((FAILED+1)); continue; }
    end2=$(date +%s)
    phase2=$((end2 - start2))
    total=$((phase1 + phase2))
    echo "$CHROM,$BACKEND,$phase1,$phase2,$total" >> "$TIMINGS_CSV"
    echo "  Phase 1: ${phase1}s  |  Phase 2: ${phase2}s  |  Total: ${total}s"
    echo ""
    DONE=$((DONE+1))
  done
done

echo ""
echo "╔══════════════════════════════════════════════════════════╗"
echo "║  Benchmark complete                                     ║"
echo "║  Successful: $DONE / $((N_CHROMS * 2))  runs"
if [[ $FAILED -gt 0 ]]; then
echo "║  Failed:     $FAILED"
fi
echo "║  Timings:    $TIMINGS_CSV"
echo "╚══════════════════════════════════════════════════════════╝"
echo ""
cat "$TIMINGS_CSV"

if [[ "$SKIP_PLOT" == "1" ]]; then
  echo "SKIP_PLOT=1 — not generating plot."
  exit 0
fi

echo ""
echo "Generating plot..."
python3 "$SCRIPT_DIR/plot_voting_benchmark.py" "$TIMINGS_CSV" "$PLOT_PATH" || {
  echo "Plot failed (install matplotlib: pip install 'ancify[evaluate]')."
  exit 0
}
echo "Plot saved: $PLOT_PATH"
