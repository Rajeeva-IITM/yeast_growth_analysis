#!/bin/bash

# Run script — classification on the bloom_256 datasets, sum gene encoding: total
# mutational load, the per-gene effect encodings added.
# Same 6014 gene columns as max, so runtime is comparable.
#
# Run from inside a `pixi shell` (uses plain `python`). Defaults to the sigma_0.5
# threshold family; pass 1.0 as the first argument for the other one:
#
#     ./scripts/run_encodings_sum.sh 1.0

alpha="${1:-0.5}"

# Classification is already the config default (objective=binary, metric=roc_auc_score,
# n_trials=100), and conf_encodings.yaml derives every path from ${alpha} and the
# selected dataset, so there is nothing left to override here.
export extra_params="alpha=${alpha}"

python src/tune_model.py --config-name=conf_encodings --multirun dataset=bloom2013_sum \
 'seed=1,2,3,4,5' $extra_params

echo "Sum_Bloom2013 (sigma_${alpha}) Done"

python src/tune_model.py --config-name=conf_encodings --multirun dataset=bloom2015_sum \
 'seed=1,2,3,4,5' $extra_params

echo "Sum_Bloom2015 (sigma_${alpha}) Done"

python src/tune_model.py --config-name=conf_encodings --multirun dataset=bloom2019_BYxRM_sum \
 'seed=1,2,3,4,5' $extra_params

echo "Sum_Bloom2019_BYxRM (sigma_${alpha}) Done"

# python src/tune_model.py --config-name=conf_encodings --multirun dataset=bloom2019_BYxM22_sum \
#  'seed=1,2,3,4,5' $extra_params
#
# echo "Sum_Bloom2019_BYxM22 (sigma_${alpha}) Done"

# python src/tune_model.py --config-name=conf_encodings --multirun dataset=bloom2019_RMxYPS163_sum \
#  'seed=1,2,3,4,5' $extra_params
#
# echo "Sum_Bloom2019_RMxYPS163 (sigma_${alpha}) Done"

echo "All Done"
