#!/bin/bash

#the next few commands are required only for actinium
export PATH=$PATH:/home/rajeeva/.local/bin/
cd /data/rajeeva/Boosting-yeast_growth_pred/ || exit
# micromamba activate boost_env
export PYTHON=~/micromamba/envs/boost_env/bin/python

# python src/tune_model.py --multirun data.path=/home/rajeeva/Project/data/regression_data/new_latent/bloom2015_regression_std.feather \
#  data.savename=Full_Bloom2013 'data.savedir=${oc.env:RUN_DIR}/regression/${data.savename}_${seed}_${run_type}' \
#  seed=1,2,3,4,5 model_params.objective=mse metric._target_=sklearn.metrics.r2_score n_trials=50 run_type=geno_only

$PYTHON src/tune_model.py --multirun 'data.path=${oc.env:DATA_DIR}/regression_data/new_latent/bloom2015_regression_std.feather' \
 data.savename=Full_Bloom2015 'data.savedir=${oc.env:RUN_DIR}/regression/${data.savename}_${seed}_full' \
 seed=1,2,3,4,5 model_params.objective=mse metric._target_=sklearn.metrics.r2_score n_trials=50
