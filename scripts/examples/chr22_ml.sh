#!/usr/bin/env bash
# Chr22, LightGBM classifier. Train a model first:
#   ancify train -c config.yaml -o model.lgb
# Then: ML_MODEL_PATH=model.lgb ./scripts/examples/chr22_ml.sh
set -euo pipefail
export METHOD=ml
if [[ -z "${ML_MODEL_PATH:-}" ]]; then
  echo "Usage: ML_MODEL_PATH=path/to/model.lgb $0"
  exit 1
fi
exec "$(dirname "$0")/../run_hg38.sh"
