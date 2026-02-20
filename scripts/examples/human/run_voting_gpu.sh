#!/usr/bin/env bash
# Human chr22 — voting, GPU backend.
set -euo pipefail
CHROM=chr22 BACKEND=gpu exec "$(dirname "$0")/../../run_hg38.sh"
