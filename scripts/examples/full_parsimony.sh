#!/usr/bin/env bash
# Full genome, Fitch parsimony on (((bonobo,chimp),gorilla),macaque).
set -euo pipefail
export CHROM=all
export METHOD=parsimony
export ANCIFY_CPUS="${ANCIFY_CPUS:-8}"
exec "$(dirname "$0")/../run_hg38.sh"
