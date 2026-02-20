#!/usr/bin/env bash
# Chr22, two-tier voting, GPU backend. Requires PyTorch + CUDA.
set -euo pipefail
export BACKEND=gpu
exec "$(dirname "$0")/../run_hg38.sh"
