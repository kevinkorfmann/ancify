#!/usr/bin/env bash
# ──────────────────────────────────────────────────────────────
# Run voting, parsimony, and likelihood on a single chromosome
# and compare the output. Results go into subdirs for easy diff.
#
# Env:
#   CHROM       chromosome to test (default: chr22)
#   BACKEND     auto | cpu | gpu (default: auto)
#   ANCIFY_CPUS workers (default: 4)
#   WORK_DIR    working directory (default: repo/ancify_compare)
# ──────────────────────────────────────────────────────────────
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

export CHROM="${CHROM:-chr22}"
export BACKEND="${BACKEND:-auto}"
export ANCIFY_CPUS="${ANCIFY_CPUS:-4}"
export WORK_DIR="${WORK_DIR:-$REPO_DIR/ancify_compare}"

METHODS=(voting parsimony likelihood)

echo "╔═══════════════════════════════════════════╗"
echo "║  Compare all methods on $CHROM"
echo "║  Backend: $BACKEND   CPUs: $ANCIFY_CPUS"
echo "║  Output:  $WORK_DIR/output_${CHROM}/"
echo "╚═══════════════════════════════════════════╝"
echo ""

for method in "${METHODS[@]}"; do
  echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
  echo "  Method: $method"
  echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
  export METHOD="$method"
  "$SCRIPT_DIR/run_hg38.sh"
  echo ""
done

# ── Compare output ───────────────────────────────────────────
OUT_BASE="$WORK_DIR/output_${CHROM}"

echo "╔═══════════════════════════════════════════╗"
echo "║  Comparison summary: $CHROM"
echo "╚═══════════════════════════════════════════╝"
echo ""
echo "Output structure:"
echo "  $OUT_BASE/"
for method in "${METHODS[@]}"; do
  fa="$OUT_BASE/$method/${CHROM}.fa"
  if [[ -f "$fa" ]]; then
    size=$(wc -c < "$fa")
    echo "  ├── $method/${CHROM}.fa  ($size bytes)"
  else
    echo "  ├── $method/${CHROM}.fa  (MISSING)"
  fi
done
echo ""

# Quick stats via Python (if ancify importable)
python3 - "$OUT_BASE" "$CHROM" "${METHODS[@]}" <<'PYEOF'
import sys
from pathlib import Path

out_base, chrom = sys.argv[1], sys.argv[2]
methods = sys.argv[3:]

seqs = {}
for m in methods:
    fa = Path(out_base) / m / f"{chrom}.fa"
    if not fa.exists():
        print(f"  {m}: output missing, skipping")
        continue
    lines = fa.read_text().splitlines()
    seqs[m] = "".join(l for l in lines if not l.startswith(">"))

if len(seqs) < 2:
    print("Need at least 2 methods to compare.")
    sys.exit(0)

for m, seq in seqs.items():
    total = len(seq)
    high = sum(1 for c in seq if c in "ACGT")
    low  = sum(1 for c in seq if c in "acgt")
    unresolved = sum(1 for c in seq if c == "n")
    missing = sum(1 for c in seq if c == "N")
    print(f"  {m:12s}  high={high/total:.1%}  low={low/total:.1%}  "
          f"unresolved={unresolved/total:.1%}  missing={missing/total:.1%}")

names = list(seqs.keys())
print()
for i in range(len(names)):
    for j in range(i + 1, len(names)):
        a, b = seqs[names[i]], seqs[names[j]]
        length = min(len(a), len(b))
        agree = sum(1 for k in range(length)
                    if a[k].upper() == b[k].upper() and a[k].upper() in "ACGTN")
        both_called = sum(1 for k in range(length)
                         if a[k].upper() in "ACGT" and b[k].upper() in "ACGT")
        if both_called > 0:
            rate = sum(1 for k in range(length)
                       if a[k].upper() in "ACGT" and b[k].upper() in "ACGT"
                       and a[k].upper() == b[k].upper()) / both_called
            print(f"  {names[i]} vs {names[j]}:  "
                  f"agreement={rate:.4%} (over {both_called:,} shared sites)")
        else:
            print(f"  {names[i]} vs {names[j]}:  no shared called sites")
PYEOF

# ── Visualize ────────────────────────────────────────────────
echo ""
echo "Generating figures..."
python3 "$SCRIPT_DIR/visualize_comparison.py" "$OUT_BASE" "$CHROM" "${METHODS[@]}" || {
  echo "(visualization skipped — install matplotlib: pip install 'ancify[evaluate]')"
}

echo ""
echo "Done. FASTA files in $OUT_BASE/{voting,parsimony,likelihood}/"
echo "      Figures in $OUT_BASE/figures/"
