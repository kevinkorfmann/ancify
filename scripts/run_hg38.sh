#!/usr/bin/env bash
# ──────────────────────────────────────────────────────────────
# Unified ancify test runner for human (hg38).
# Downloads UCSC data and runs the pipeline with configurable
# scope, method, and backend.
#
# Env:
#   CHROM       chr22 (default) or "all" (autosomes + chrX)
#   METHOD      voting (default), parsimony, likelihood, ml
#   BACKEND     auto (default), cpu, gpu
#   ANCIFY_CPUS number of workers (default: 4)
#   WORK_DIR    working directory (default: repo/ancify_test)
#   ML_MODEL_PATH  path to trained .lgb model (required if ml)
# ──────────────────────────────────────────────────────────────
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

CHROM="${CHROM:-chr22}"
METHOD="${METHOD:-voting}"
BACKEND="${BACKEND:-auto}"
WORK_DIR="${WORK_DIR:-$REPO_DIR/ancify_test}"
BASE_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38"

echo "═══════════════════════════════════════════"
echo "  ancify test run"
echo "  scope:   $CHROM"
echo "  method:  $METHOD"
echo "  backend: $BACKEND"
echo "  workdir: $WORK_DIR"
echo "═══════════════════════════════════════════"

mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

# ── Download chromosome lengths ──────────────────────────────
if [[ ! -f hg38.chrom.sizes ]]; then
  echo "[download] hg38 chromosome sizes..."
  curl -sL "$BASE_URL/bigZips/hg38.chrom.sizes" -o hg38.chrom.sizes
fi

if [[ "$CHROM" == all ]]; then
  cp hg38.chrom.sizes chromoLens.txt
else
  awk -v c="$CHROM" '$1 == c' hg38.chrom.sizes > chromoLens.txt
  if [[ ! -s chromoLens.txt ]]; then
    echo "ERROR: chromosome '$CHROM' not found in hg38.chrom.sizes"
    exit 1
  fi
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

if [[ "$CHROM" == all ]]; then
  download hg38.panPan3.net.axt.gz   vsPanPan3
  download hg38.panTro6.net.axt.gz   vsPanTro6
  download hg38.gorGor6.net.axt.gz   vsGorGor6
  download hg38.rheMac10.net.axt.gz  vsRheMac10
  INNER_NAMES="bonobo chimp gorilla"
  OUTER_NAME="macaque"
  TREE_PARSIMONY='(((bonobo,chimp),gorilla),macaque)'
  TREE_LIKELIHOOD='(((bonobo:0.008,chimp:0.008):0.002,gorilla:0.009):0.020,macaque:0.038)'
else
  download hg38.panTro6.net.axt.gz   vsPanTro6
  download hg38.rheMac10.net.axt.gz  vsRheMac10
  INNER_NAMES="chimp"
  OUTER_NAME="macaque"
  TREE_PARSIMONY='(chimp,macaque)'
  TREE_LIKELIHOOD='(chimp:0.008,macaque:0.038)'
fi

# ── Validate ML ──────────────────────────────────────────────
if [[ "$METHOD" == ml && -z "${ML_MODEL_PATH:-}" ]]; then
  echo "ERROR: METHOD=ml requires ML_MODEL_PATH=/path/to/model.lgb"
  echo "Train first: ancify train -c config.yaml -o model.lgb"
  exit 1
fi

# ── Build output dir name ────────────────────────────────────
if [[ "$METHOD" == voting ]]; then
  OUT_DIR="$WORK_DIR/output_${CHROM}/voting"
else
  OUT_DIR="$WORK_DIR/output_${CHROM}/${METHOD}"
fi

# ── Write config ─────────────────────────────────────────────
CONFIG="$WORK_DIR/config_${CHROM}_${METHOD}.yaml"

inner_block() {
  for name in $INNER_NAMES; do
    local file
    case "$name" in
      bonobo)  file=hg38.panPan3.net.axt.gz ;;
      chimp)   file=hg38.panTro6.net.axt.gz ;;
      gorilla) file=hg38.gorGor6.net.axt.gz ;;
    esac
    echo "    - name: $name"
    echo "      alignment: $WORK_DIR/$file"
  done
}

outer_file() {
  case "$OUTER_NAME" in
    macaque) echo "hg38.rheMac10.net.axt.gz" ;;
  esac
}

{
  cat <<EOF
focal_species: human
chromosome_lengths: $WORK_DIR/chromoLens.txt
EOF

  if [[ "$CHROM" != all ]]; then
    echo "chromosomes:"
    echo "  - $CHROM"
  else
    echo "chromosomes:"
    for c in chr{1..22} chrX; do echo "  - $c"; done
  fi

  echo "outgroups:"
  echo "  inner:"
  inner_block
  echo "  outer:"
  echo "    - name: $OUTER_NAME"
  echo "      alignment: $WORK_DIR/$(outer_file)"

  cat <<EOF
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
