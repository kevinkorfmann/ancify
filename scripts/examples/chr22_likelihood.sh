#!/usr/bin/env bash
# Chr22, likelihood (Felsenstein pruning, HKY85 model).
set -euo pipefail
export METHOD=likelihood
exec "$(dirname "$0")/../run_hg38.sh"
