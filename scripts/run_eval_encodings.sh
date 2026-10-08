#!/bin/bash

# Predictions + per-condition classification metrics for the bloom_256 encoding runs.
# Run from inside a `pixi shell`, from the project root. Defaults to max at sigma_0.5:
#
#     ./scripts/run_eval_encodings.sh              # max, sigma_0.5
#     ./scripts/run_eval_encodings.sh sum 0.5      # sum
#
# compare_performance.py loops models x datasets, so one invocation writes 9 prediction
# parquets and 9 metric csvs per encoding: the 3 matched pairs plus the 6 cross-dataset
# combinations (a Bloom2013-trained model scored on Bloom2015 data, and so on).

set -euo pipefail

encoding="${1:-max}"

case "${encoding}" in
  max)   prefix=Max ;;
  sum)   prefix=Sum ;;
  count) prefix=Count ;;
  *) echo "unknown encoding '${encoding}' (expected max, sum or count)" >&2; exit 1 ;;
esac

alpha="${2:-0.5}"

# Which family of training runs to score. No run in this project uses linear trees any
# more; the name is kept because it is baked into every output path.
variant="${3:-no_linear_tree}"

# PROJECT_DIR lives in .env, which only python-dotenv reads -- the shell needs it too, to
# create out_path below. compare_performance.py writes straight into it and does not mkdir.
set -a; . ./.env; set +a

out_path="${PROJECT_DIR}/performance/classification_encodings/${variant}/sigma_${alpha}/${encoding}"
mkdir -p "${out_path}"

export extra_params="encoding=${encoding} encoding_prefix=${prefix} alpha=${alpha} variant=${variant}"

python src/compare_performance.py --config-name=eval_encodings $extra_params

echo "Eval ${prefix} (sigma_${alpha}) Done -> ${out_path}"
