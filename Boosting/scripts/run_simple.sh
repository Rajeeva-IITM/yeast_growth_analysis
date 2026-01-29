#!/bin/bash

# export python=/data/rajeeva/micromamba/envs/boost_env/bin/python
# cd "/data/rajeeva/Boosting-yeast_growth_pred/" ||
# export seed=1
export SEED="1,2,3,4,5"   # Suffix of the directory created
export extra_params='data.savedir=${oc.env:ROOT_DIR}/runs/classification_other_models/${data.savename}_${seed}_${run_type}
           n_trials=50'

python src/tune_and_train_simple_models_2.py --multirun 'data.path=${oc.env:DATA_DIR}/full/varying_sigma/sigma_0.5/bloom2013_clf.feather' \
 data.savename=Full_Bloom2013 'seed=1,2,3,4,5' $extra_params

echo "Full_Bloom2013 Done"


python src/tune_and_train_simple_models_2.py --multirun 'data.path=${oc.env:DATA_DIR}/full/varying_sigma/sigma_0.5/bloom2019_clf.feather' \
 data.savename=Full_Bloom2019_BYxRM 'seed=1,2,3,4,5' $extra_params

echo "Full_Bloom2019 Done"

# python src/tune_and_train_simple_models_2.py --multirun 'data.path=${oc.env:DATA_DIR}/full/varying_sigma/sigma_0.5/bloom2013_clf.feather' \
#  data.savename=Full_Bloom2019_BYxM22 'seed=1,2,3,4,5' $extra_params

# echo "Full_Bloom2019_BYxM22 Done"

# python src/tune_and_train_simple_models_2.py --multirun 'data.path=${oc.env:DATA_DIR}/full/varying_sigma/sigma_0.5/bloom2013_clf.feather' \
#  data.savename=Full_Bloom2019_RMxYPS163 'seed=1,2,3,4,5' $extra_params

# echo "Full_Bloom2019_RMxYPS163 Done"

python src/tune_and_train_simple_models_2.py --multirun 'data.path=${oc.env:DATA_DIR}/full/varying_sigma/sigma_0.5/bloom2015_clf.feather' \
 data.savename=Full_Bloom2015 'seed=1,2,3,4,5' $extra_params

echo "Full_Bloom2015 Done"

echo "All Done"
