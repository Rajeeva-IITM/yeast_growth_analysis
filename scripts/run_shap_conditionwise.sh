#!/bin/bash

# Per-condition SHAP over the bloom_256 encoding runs. Run from inside a `pixi shell`, from
# the project root. Defaults to the max encoding at sigma_0.5, exact method:
#
#     ./scripts/run_shap_conditionwise.sh                              # max, sigma_0.5, exact
#     ./scripts/run_shap_conditionwise.sh sum 0.5                      # sum
#     ./scripts/run_shap_conditionwise.sh max 0.5 no_linear_tree saabas  # the fast companion
#
# One invocation covers all three datasets and every condition in each: 39 for Bloom2013, 18
# for Bloom2015, 31 for Bloom2019_BYxRM, each explained under all five seeds.
#
# The two methods write to separate trees and must never be mixed:
#
#   exact   results/shap_conditionwise/...         exact interventional TreeSHAP, ~2.5-3 h
#   saabas  results/shap_conditionwise_saabas/...  _cext.dense_tree_saabas, a few minutes
#
# Only `exact` actually zeroes the latent block, and it does so exactly (0.000e+00 on every
# seed measured). Saabas was expected to zero it too -- a within-condition background rewrites
# every internal node expectation, which ought to make a latent split contribute nothing --
# but measurement says otherwise: 4e-2 to 9e-2 leaks through on four of five Bloom2013 seeds,
# and enlarging the background makes it worse rather than better. Report from `exact`. The
# companion exists to quantify how far the cheap estimator drifts, not to stand in for it.
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

# Which family of training runs to read. No run in this project uses linear trees any more;
# the name is kept because it is baked into every output path.
variant="${3:-no_linear_tree}"

method="${4:-exact}"

# Saabas leaks 4e-2 to 9e-2 into the latent block on four of five seeds, so it cannot meet
# the tolerance the exact run is held to; the ceiling is lifted for it and each condition's
# leak goes to the run log instead.
case "${method}" in
  exact)  approximate=false; tolerance=1e-10; tree=shap_conditionwise ;;
  saabas) approximate=true;  tolerance=1.0;   tree=shap_conditionwise_saabas ;;
  *) echo "unknown method '${method}' (expected exact or saabas)" >&2; exit 1 ;;
esac

# PROJECT_DIR lives in .env, which only python-dotenv reads -- the shell needs it too, to
# create out_path below.
set -a; . ./.env; set +a

out_path="${PROJECT_DIR}/results/${tree}/${variant}/sigma_${alpha}/${encoding}"
mkdir -p "${out_path}"

export extra_params="encoding=${encoding} encoding_prefix=${prefix} alpha=${alpha} variant=${variant} approximate=${approximate} latent_tolerance=${tolerance} out_path=${out_path}"

python src/shap_conditionwise.py --config-name=interpret_conditionwise $extra_params

echo "Conditionwise SHAP ${prefix} (sigma_${alpha}, ${method}) Done -> ${out_path}"
