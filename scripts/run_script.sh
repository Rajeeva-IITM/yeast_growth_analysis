#!/bin/bash

# Run script — classification pipeline on the sigma_0.5 datasets.
# Run from inside a `pixi shell` (uses plain `python`).

# Classification is already the config default (objective=binary, metric=roc_auc_score,
# n_trials=100), so the only override we need is the run output directory.
export extra_params='data.savedir=${oc.env:RUN_DIR}/classification_sigma/0.5_sigma/${data.savename}_${seed}_${run_type}'

python src/tune_model.py --multirun 'data.path=${oc.env:DATA_DIR}/full/varying_sigma/sigma_0.5/bloom2013_clf.feather' \
 data.savename=Full_Bloom2013 'seed=1,2,3,4,5' $extra_params

echo "Full_Bloom2013 Done"

python src/tune_model.py --multirun 'data.path=${oc.env:DATA_DIR}/full/varying_sigma/sigma_0.5/bloom2015_clf.feather' \
 data.savename=Full_Bloom2015 'seed=1,2,3,4,5' $extra_params

echo "Full_Bloom2015 Done"

python src/tune_model.py --multirun 'data.path=${oc.env:DATA_DIR}/full/varying_sigma/sigma_0.5/bloom2019_clf.feather' \
 data.savename=Full_Bloom2019_BYxRM 'seed=1,2,3,4,5' $extra_params

echo "Full_Bloom2019_BYxRM Done"

# python src/tune_model.py --multirun 'data.path=${oc.env:DATA_DIR}/full/varying_sigma/sigma_0.5/bloom2019_BYxM22_clf.feather' \
#  data.savename=Full_Bloom2019_BYxM22 'seed=1,2,3,4,5' $extra_params
#
# echo "Full_Bloom2019_BYxM22 Done"
#
# python src/tune_model.py --multirun 'data.path=${oc.env:DATA_DIR}/full/varying_sigma/sigma_0.5/bloom2019_RMxYPS163_clf.feather' \
#  data.savename=Full_Bloom2019_RMxYPS163 'seed=1,2,3,4,5' $extra_params
#
# echo "Full_Bloom2019_RMxYPS163 Done"

echo "All Done"
