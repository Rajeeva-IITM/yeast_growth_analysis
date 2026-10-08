#!/bin/bash

# Run script — classification on the bloom_256 datasets, max gene encoding: the
# largest effect encoding among a gene's variants.
# This is the encoding every model before August 2026 used, rebuilt on the 256-dim
# latent, so it is the baseline the sum and count families are read against.
#
# Run from inside a `pixi shell` (uses plain `python`). Defaults to the sigma_0.5
# threshold family; pass 1.0 as the first argument for the other one:
#
#     ./scripts/run_encodings_max.sh 1.0

alpha="${1:-0.5}"

# Classification is already the config default (objective=binary, metric=roc_auc_score,
# n_trials=100), and conf_encodings.yaml derives every path from ${alpha} and the
# selected dataset, so there is nothing left to override here.
export extra_params="alpha=${alpha}"

python src/tune_model.py --config-name=conf_encodings --multirun dataset=bloom2013_max \
 'seed=1,2,3,4,5' $extra_params

echo "Max_Bloom2013 (sigma_${alpha}) Done"

python src/tune_model.py --config-name=conf_encodings --multirun dataset=bloom2015_max \
 'seed=1,2,3,4,5' $extra_params

echo "Max_Bloom2015 (sigma_${alpha}) Done"

python src/tune_model.py --config-name=conf_encodings --multirun dataset=bloom2019_BYxRM_max \
 'seed=1,2,3,4,5' $extra_params

echo "Max_Bloom2019_BYxRM (sigma_${alpha}) Done"

# python src/tune_model.py --config-name=conf_encodings --multirun dataset=bloom2019_BYxM22_max \
#  'seed=1,2,3,4,5' $extra_params
#
# echo "Max_Bloom2019_BYxM22 (sigma_${alpha}) Done"

# python src/tune_model.py --config-name=conf_encodings --multirun dataset=bloom2019_RMxYPS163_max \
#  'seed=1,2,3,4,5' $extra_params
#
# echo "Max_Bloom2019_RMxYPS163 (sigma_${alpha}) Done"

echo "All Done"
