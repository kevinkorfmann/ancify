#!/usr/bin/env bash
# ──────────────────────────────────────────────────────────────
# ancify runner for Brassica rapa (braRap1). Single chromosome.
#
# Inner outgroup: B. oleracea (same genus)
# Outer outgroup: Arabidopsis thaliana (~20 Mya, same family)
#
# NOTE: Unlike the human / mouse / Drosophila runners, this script
# does NOT auto-download alignments from UCSC.  You must provide:
#   1. braRap1.chromLens.txt      — chromosome name + length
#   2. braRap1.braOleracea.net.axt.gz  — B. rapa vs B. oleracea
#   3. braRap1.araTha1.net.axt.gz      — B. rapa vs A. thaliana
# Place them in WORK_DIR (or symlink).
#
# Env:
#   CHROM       A01 (default) — any single chromosome in the lengths file
#   METHOD      voting (default), parsimony, likelihood, ml
#   BACKEND     auto (default), cpu, gpu
#   ANCIFY_CPUS number of workers (default: 4)
#   WORK_DIR    working directory (default: repo/ancify_test_brassica)
#   ML_MODEL_PATH  path to trained .lgb model (required if ml)
# ──────────────────────────────────────────────────────────────
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

CHROM="${CHROM:-A01}"
METHOD="${METHOD:-voting}"
BACKEND="${BACKEND:-auto}"
WORK_DIR="${WORK_DIR:-$REPO_DIR/ancify_test_brassica}"

echo "═══════════════════════════════════════════"
echo "  ancify — Brassica rapa (braRap1)"
echo "  chromosome: $CHROM"
echo "  method:     $METHOD"
echo "  backend:    $BACKEND"
echo "  workdir:    $WORK_DIR"
echo "═══════════════════════════════════════════"

mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

# ── Verify required files ────────────────────────────────────
CHROM_LENS="$WORK_DIR/braRap1.chromLens.txt"
INNER_ALN="$WORK_DIR/braRap1.braOleracea.net.axt.gz"
OUTER_ALN="$WORK_DIR/braRap1.araTha1.net.axt.gz"

for f in "$CHROM_LENS" "$INNER_ALN" "$OUTER_ALN"; do
  if [[ ! -f "$f" ]]; then
    echo "ERROR: required file not found: $f"
    echo "Place alignment data in WORK_DIR before running. See header comments."
    exit 1
  fi
done

awk -v c="$CHROM" '$1 == c' "$CHROM_LENS" > chromoLens.txt
if [[ ! -s chromoLens.txt ]]; then
  echo "ERROR: chromosome '$CHROM' not found in $CHROM_LENS"
  exit 1
fi

TREE_PARSIMONY='(brassica_oleracea,arabidopsis_thaliana)'
TREE_LIKELIHOOD='(brassica_oleracea:0.030,arabidopsis_thaliana:0.200)'

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
focal_species: brassica_rapa
chromosome_lengths: $WORK_DIR/chromoLens.txt
chromosomes:
  - $CHROM
outgroups:
  inner:
    - name: brassica_oleracea
      alignment: $INNER_ALN
  outer:
    - name: arabidopsis_thaliana
      alignment: $OUTER_ALN
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
