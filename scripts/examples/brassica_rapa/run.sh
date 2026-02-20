#!/usr/bin/env bash
# Brassica rapa chromosome A01 — voting (default method), CPU.
#
# Requires alignment data in WORK_DIR before running — see scripts/run_brassica.sh
# for the list of required files.
set -euo pipefail
CHROM=A01 exec "$(dirname "$0")/../../run_brassica.sh"
