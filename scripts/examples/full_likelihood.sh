#!/usr/bin/env bash
# Full genome, likelihood (Felsenstein pruning, HKY85).
set -euo pipefail
export CHROM=all
export METHOD=likelihood
export ANCIFY_CPUS="${ANCIFY_CPUS:-8}"
exec "$(dirname "$0")/../run_hg38.sh"
