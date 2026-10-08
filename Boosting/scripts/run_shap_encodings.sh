#!/bin/bash

# SHAP analysis over the bloom_256 encoding runs. Run from inside a `pixi shell`, from the
# project root. Defaults to the max encoding at sigma_0.5:
#
#     ./scripts/run_shap_encodings.sh              # max, sigma_0.5
#     ./scripts/run_shap_encodings.sh sum 0.5      # sum
#     ./scripts/run_shap_encodings.sh max 1.0      # max at the other threshold
#
# One invocation covers all three datasets: shap_analysis.py loops over
# model_load_keys.model_names internally, and for each one globs all five seeds and
# averages the per-feature mean |SHAP| across them.
#
# Needs `shap` in the environment (shap 0.45.1 verified).

set -euo pipefail

encoding="${1:-max}"

case "${encoding}" in
  max)   prefix=Max ;;
  sum)   prefix=Sum ;;
  count) prefix=Count ;;
  *) echo "unknown encoding '${encoding}' (expected max, sum or count)" >&2; exit 1 ;;
esac

alpha="${2:-0.5}"

# Which family of training runs to read. No run in this project uses linear trees any
# more; the name is kept because it is baked into every output path. The superseded
# linear-tree results are still on disk one level shallower, but their models were deleted
# and cannot be regenerated, so nothing here should point at them.
variant="${3:-no_linear_tree}"


# PROJECT_DIR lives in .env, which only python-dotenv reads -- the shell needs it too, to
# create out_path below. shap_analysis.py writes straight into out_path and does not mkdir it.
set -a; . ./.env; set +a

out_path="${PROJECT_DIR}/results/shap_encodings/${variant}/sigma_${alpha}/${encoding}"
mkdir -p "${out_path}"

export extra_params="encoding=${encoding} encoding_prefix=${prefix} alpha=${alpha} variant=${variant}"

python src/shap_analysis.py --config-name=interpret_encodings $extra_params

echo "SHAP ${prefix} (sigma_${alpha}) Done -> ${out_path}"
