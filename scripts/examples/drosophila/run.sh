#!/usr/bin/env bash
# Drosophila melanogaster chromosome 4 — voting (default method), CPU.
# Chrom 4 is the smallest (~1.3 Mb) — fast to download and run.
set -euo pipefail
CHROM=4 exec "$(dirname "$0")/../../run_dm6.sh"
