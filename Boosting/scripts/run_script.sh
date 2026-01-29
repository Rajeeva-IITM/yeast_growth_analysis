#!/bin/bash

# cd "/home/rajeeva/Project/boosting/" ||
# Run script

export python=/data/rajeeva/micromamba/envs/boost_env/bin/python
cd "/data/rajeeva/Boosting-yeast_growth_pred/" ||
# export seed=1
export SEED="1,2,3,4,5"   # Suffix of the directory created
export extra_params='data.savedir=${oc.env:RUN_DIR}/regression/minmax/${data.savename}_${seed}_${run_type} n_trials=50 model_params.objective=mse
           metric._target_=sklearn.metrics.r2_score'

$python src/tune_model.py --multirun 'data.path=${oc.env:DATA_DIR}/regression_data/minmax_scaled/bloom2013_regression_minmax.feather' \
 data.savename=Full_Bloom2013 'seed=1,2,3,4,5' $extra_params

echo "Full_Bloom2013 Done"

$python src/tune_model.py --multirun 'data.path=${oc.env:DATA_DIR}/regression_data/minmax_scaled/bloom2015_regression_minmax.feather' \
 data.savename=Full_Bloom2015 'seed=1,2,3,4,5' $extra_params

echo "Full_Bloom2015 Done"

$python src/tune_model.py --multirun 'data.path=${oc.env:DATA_DIR}/regression_data/minmax_scaled/bloom2019_regression_minmax.feather' \
 data.savename=Full_Bloom2019_BYxRM 'seed=2,3,4,5' $extra_params

echo "Full_Bloom2019 Done"

$python src/tune_model.py --multirun 'data.path=${oc.env:DATA_DIR}/regression_data/minmax_scaled/bloom2019_BYxM22_regression_minmax.feather' \
 data.savename=Full_Bloom2019_BYxM22 'seed=1,2,3,4,5' $extra_params

echo "Full_Bloom2019_BYxM22 Done"

$python src/tune_model.py --multirun 'data.path=${oc.env:DATA_DIR}/regression_data/minmax_scaled/bloom2019_RMxYPS163_regression_minmax.feather' \
 data.savename=Full_Bloom2019_RMxYPS163 'seed=1,2,3,4,5' $extra_params

echo "Full_Bloom2019_RMxYPS163 Done"

echo "All Done"
