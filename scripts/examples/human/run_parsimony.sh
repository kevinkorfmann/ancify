#!/usr/bin/env bash
# Human chr22 — Fitch parsimony.
set -euo pipefail
CHROM=chr22 METHOD=parsimony exec "$(dirname "$0")/../../run_hg38.sh"
