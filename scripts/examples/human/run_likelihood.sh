#!/usr/bin/env bash
# Human chr22 — likelihood (HKY85).
set -euo pipefail
CHROM=chr22 METHOD=likelihood exec "$(dirname "$0")/../../run_hg38.sh"
