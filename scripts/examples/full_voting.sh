#!/usr/bin/env bash
# Full genome (chr1-22 + chrX), voting, 4 outgroups.
set -euo pipefail
export CHROM=all
export ANCIFY_CPUS="${ANCIFY_CPUS:-8}"
exec "$(dirname "$0")/../run_hg38.sh"
