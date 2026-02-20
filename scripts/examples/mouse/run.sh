#!/usr/bin/env bash
# Mouse chr19 — voting (default method), CPU.
set -euo pipefail
CHROM=chr19 exec "$(dirname "$0")/../../run_mm39.sh"
