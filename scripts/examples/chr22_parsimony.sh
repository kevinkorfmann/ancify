#!/usr/bin/env bash
# Chr22, Fitch parsimony on (chimp, macaque) tree.
set -euo pipefail
export METHOD=parsimony
exec "$(dirname "$0")/../run_hg38.sh"
