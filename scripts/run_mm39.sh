#!/usr/bin/env bash
# ──────────────────────────────────────────────────────────────
# ancify runner for mouse (mm39). Single chromosome.
# Downloads UCSC data and runs the pipeline.
#
# Inner outgroup: rat (rn7, ~12 Mya)
# Outer outgroup: rabbit (oryCun2, ~90 Mya)
#
# Env:
#   CHROM       chr19 (default) — any single mm39 chromosome
#   METHOD      voting (default), parsimony, likelihood, ml
#   BACKEND     auto (default), cpu, gpu
#   ANCIFY_CPUS number of workers (default: 4)
#   WORK_DIR    working directory (default: repo/ancify_test_mouse)
#   ML_MODEL_PATH  path to trained .lgb model (required if ml)
# ──────────────────────────────────────────────────────────────
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

CHROM="${CHROM:-chr19}"
METHOD="${METHOD:-voting}"
BACKEND="${BACKEND:-auto}"
WORK_DIR="${WORK_DIR:-$REPO_DIR/ancify_test_mouse}"
BASE_URL="https://hgdownload.soe.ucsc.edu/goldenPath/mm39"

echo "═══════════════════════════════════════════"
echo "  ancify — mouse (mm39)"
echo "  chromosome: $CHROM"
echo "  method:     $METHOD"
echo "  backend:    $BACKEND"
echo "  workdir:    $WORK_DIR"
echo "═══════════════════════════════════════════"

mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

# ── Download chromosome lengths ──────────────────────────────
if [[ ! -f mm39.chrom.sizes ]]; then
  echo "[download] mm39 chromosome sizes..."
  curl -sL "$BASE_URL/bigZips/mm39.chrom.sizes" -o mm39.chrom.sizes
fi

awk -v c="$CHROM" '$1 == c' mm39.chrom.sizes > chromoLens.txt
if [[ ! -s chromoLens.txt ]]; then
  echo "ERROR: chromosome '$CHROM' not found in mm39.chrom.sizes"
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

download mm39.rn7.net.axt.gz      vsRn7
download mm39.oryCun2.net.axt.gz  vsOryCun2

TREE_PARSIMONY='(rat,rabbit)'
TREE_LIKELIHOOD='(rat:0.012,rabbit:0.090)'

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
focal_species: mouse
chromosome_lengths: $WORK_DIR/chromoLens.txt
chromosomes:
  - $CHROM
outgroups:
  inner:
    - name: rat
      alignment: $WORK_DIR/mm39.rn7.net.axt.gz
  outer:
    - name: rabbit
      alignment: $WORK_DIR/mm39.oryCun2.net.axt.gz
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
