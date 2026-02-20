#!/usr/bin/env bash
# ──────────────────────────────────────────────────────────────
# ancify runner for Drosophila melanogaster (dm6). Single chromosome.
# Downloads UCSC data and runs the pipeline.
#
# Inner outgroups: D. simulans (droSim2), D. sechellia (droSec1)
# Outer outgroup:  D. yakuba (droYak3, ~6 Mya)
#
# Env:
#   CHROM       4 (default) — any single dm6 chromosome (2L, 2R, 3L, 3R, 4, X)
#   METHOD      voting (default), parsimony, likelihood, ml
#   BACKEND     auto (default), cpu, gpu
#   ANCIFY_CPUS number of workers (default: 4)
#   WORK_DIR    working directory (default: repo/ancify_test_dm6)
#   ML_MODEL_PATH  path to trained .lgb model (required if ml)
# ──────────────────────────────────────────────────────────────
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

CHROM="${CHROM:-4}"
METHOD="${METHOD:-voting}"
BACKEND="${BACKEND:-auto}"
WORK_DIR="${WORK_DIR:-$REPO_DIR/ancify_test_dm6}"
BASE_URL="https://hgdownload.soe.ucsc.edu/goldenPath/dm6"

echo "═══════════════════════════════════════════"
echo "  ancify — Drosophila melanogaster (dm6)"
echo "  chromosome: $CHROM"
echo "  method:     $METHOD"
echo "  backend:    $BACKEND"
echo "  workdir:    $WORK_DIR"
echo "═══════════════════════════════════════════"

mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

# ── Download chromosome lengths ──────────────────────────────
if [[ ! -f dm6.chrom.sizes ]]; then
  echo "[download] dm6 chromosome sizes..."
  curl -sL "$BASE_URL/bigZips/dm6.chrom.sizes" -o dm6.chrom.sizes
fi

awk -v c="$CHROM" '$1 == c' dm6.chrom.sizes > chromoLens.txt
if [[ ! -s chromoLens.txt ]]; then
  echo "ERROR: chromosome '$CHROM' not found in dm6.chrom.sizes"
  exit 1
fi

# ── Download alignments ──────────────────────────────────────
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

download dm6.droSim2.net.axt.gz  vsDroSim2
download dm6.droSec1.net.axt.gz  vsDroSec1
download dm6.droYak3.net.axt.gz  vsDroYak3

TREE_PARSIMONY='((simulans,sechellia),yakuba)'
TREE_LIKELIHOOD='((simulans:0.004,sechellia:0.004):0.002,yakuba:0.006)'

# ── Validate ML ──────────────────────────────────────────────
if [[ "$METHOD" == ml && -z "${ML_MODEL_PATH:-}" ]]; then
  echo "ERROR: METHOD=ml requires ML_MODEL_PATH=/path/to/model.lgb"
  exit 1
fi

OUT_DIR="$WORK_DIR/output_${CHROM}/${METHOD}"

# ── Write config ─────────────────────────────────────────────
CONFIG="$WORK_DIR/config_${CHROM}_${METHOD}.yaml"

{
  cat <<EOF
focal_species: drosophila_melanogaster
chromosome_lengths: $WORK_DIR/chromoLens.txt
chromosomes:
  - $CHROM
outgroups:
  inner:
    - name: simulans
      alignment: $WORK_DIR/dm6.droSim2.net.axt.gz
    - name: sechellia
      alignment: $WORK_DIR/dm6.droSec1.net.axt.gz
  outer:
    - name: yakuba
      alignment: $WORK_DIR/dm6.droYak3.net.axt.gz
work_dir: $WORK_DIR
output_dir: $OUT_DIR
num_cpus: ${ANCIFY_CPUS:-4}
method: $METHOD
backend: $BACKEND
EOF

  case "$METHOD" in
    parsimony)
      echo "tree: \"$TREE_PARSIMONY\""
      ;;
    likelihood)
      echo "tree: \"$TREE_LIKELIHOOD\""
      echo "substitution_model: HKY85"
      echo "model_kappa: 2.0"
      ;;
    ml)
      echo "ml_model_path: $ML_MODEL_PATH"
      ;;
  esac
} > "$CONFIG"

echo "[config] $CONFIG"

# ── Run ──────────────────────────────────────────────────────
echo "[run] ancify run -c $CONFIG"
cd "$REPO_DIR"
ancify run -c "$CONFIG"

echo ""
echo "Done. Output: $OUT_DIR/"
