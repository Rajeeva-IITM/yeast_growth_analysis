#!/bin/bash

# Run script — classification on the bloom_256 datasets, count gene encoding: three
# columns per gene holding how many variants of effect
# score 1, 2 and 3 the strain carries.
#
# NOTE: 18042 gene columns, 3x every other family. These runs dominate the compute
# bill and Count_Bloom2015 peaks around 40-50 GB resident, so do not run two of them
# at once. Append n_trials=<n> to the command line to cut the Optuna budget.
#
# Run from inside a `pixi shell` (uses plain `python`). Defaults to the sigma_0.5
# threshold family; pass 1.0 as the first argument for the other one:
#
#     ./scripts/run_encodings_count.sh 1.0

alpha="${1:-0.5}"

# Classification is already the config default (objective=binary, metric=roc_auc_score,
# n_trials=100), and conf_encodings.yaml derives every path from ${alpha} and the
# selected dataset, so there is nothing left to override here.
export extra_params="alpha=${alpha}"

python src/tune_model.py --config-name=conf_encodings --multirun dataset=bloom2013_count \
 'seed=1,2,3,4,5' $extra_params

echo "Count_Bloom2013 (sigma_${alpha}) Done"

python src/tune_model.py --config-name=conf_encodings --multirun dataset=bloom2015_count \
 'seed=1,2,3,4,5' $extra_params

echo "Count_Bloom2015 (sigma_${alpha}) Done"

python src/tune_model.py --config-name=conf_encodings --multirun dataset=bloom2019_BYxRM_count \
 'seed=1,2,3,4,5' $extra_params

echo "Count_Bloom2019_BYxRM (sigma_${alpha}) Done"

# python src/tune_model.py --config-name=conf_encodings --multirun dataset=bloom2019_BYxM22_count \
#  'seed=1,2,3,4,5' $extra_params
#
# echo "Count_Bloom2019_BYxM22 (sigma_${alpha}) Done"

# python src/tune_model.py --config-name=conf_encodings --multirun dataset=bloom2019_RMxYPS163_count \
#  'seed=1,2,3,4,5' $extra_params
#
# echo "Count_Bloom2019_RMxYPS163 (sigma_${alpha}) Done"

echo "All Done"
