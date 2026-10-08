#!/bin/bash

# Partial-dependence / interaction analysis over the bloom_256 encoding runs. Run from
# inside a `pixi shell`, from the project root. Defaults to the max encoding at sigma_0.5:
#
#     ./scripts/run_interaction_encodings.sh              # max, sigma_0.5
#     ./scripts/run_interaction_encodings.sh sum 0.5      # sum
#     ./scripts/run_interaction_encodings.sh max 1.0      # max at the other threshold
#
# One invocation covers all three datasets and all five seeds each. No statistic reported
# here comes from shap -- see the header of src/interaction_analysis.py for why (SHAP is an
# observational attribution, partial dependence is interventional, and the ranking has to
# match the curves). It does read the mean|SHAP| csvs, if present, to widen the candidate
# pool; run scripts/run_shap_encodings.sh for the same encoding and alpha first, or the
# pool falls back to split gain alone.

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
# create out_path below.
set -a; . ./.env; set +a

out_path="${PROJECT_DIR}/results/interactions/${variant}/sigma_${alpha}/${encoding}/pdp"
mkdir -p "${out_path}"

export extra_params="encoding=${encoding} encoding_prefix=${prefix} alpha=${alpha} variant=${variant}"

python src/interaction_analysis.py --config-name=interaction $extra_params

echo "Interactions ${prefix} (sigma_${alpha}) Done -> ${out_path}"
