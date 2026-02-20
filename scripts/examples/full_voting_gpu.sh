#!/usr/bin/env bash
# Full genome, voting, GPU backend. Requires PyTorch + CUDA.
set -euo pipefail
export CHROM=all
export BACKEND=gpu
export ANCIFY_CPUS="${ANCIFY_CPUS:-8}"
exec "$(dirname "$0")/../run_hg38.sh"
