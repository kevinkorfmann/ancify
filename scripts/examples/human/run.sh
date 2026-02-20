#!/usr/bin/env bash
# Human chr22 — voting (default method), CPU.
set -euo pipefail
CHROM=chr22 exec "$(dirname "$0")/../../run_hg38.sh"
