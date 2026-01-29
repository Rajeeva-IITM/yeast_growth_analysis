#!/bin/bash

# cd /home/rajeeva/Project/boosting/

export train_data='${oc.env:DATA_DIR}/regression_data/new_latent/bloom2019_regression_std.feather'
export train_name=Bloom2019_BYxRM
export prefix="Full_"
export suffix="robust_samples"
export extra_params='data.savedir=${oc.env:RUN_DIR}/regression_sampling/${data.savename}/'


python src/data_sampling_analysis.py data.train_data=${train_data} \
 data.savename="${prefix}${train_name}_on_Bloom2013_${suffix}" \
 data.test_data='${oc.env:DATA_DIR}/regression_data/new_latent/bloom2013_regression_std.feather' $extra_params
echo "${train_name} on Bloom2013"

python src/data_sampling_analysis.py data.train_data=${train_data} \
 data.savename="${prefix}${train_name}_on_Bloom2015_${suffix}" \
 data.test_data='${oc.env:DATA_DIR}/regression_data/new_latent/bloom2015_regression_std.feather' $extra_params
echo "${train_name} on Bloom2013"

python src/data_sampling_analysis.py data.train_data=${train_data} \
 data.savename="${prefix}${train_name}_on_Bloom2019_BYxRM_${suffix}" \
 data.test_data='${oc.env:DATA_DIR}/regression_data/new_latent/bloom2019_regression_std.feather' $extra_params
echo "${train_name} on Bloom2019_new"

python src/data_sampling_analysis.py data.train_data=${train_data} \
 data.savename="${prefix}${train_name}_on_Bloom2019_BYxM22_${suffix}" \
 data.test_data='${oc.env:DATA_DIR}/regression_data/new_latent/bloom2019_BYxM22_regression_std.feather' $extra_params
echo "${train_name} on Bloom2019_BYxM22"

python src/data_sampling_analysis.py data.train_data=${train_data} \
 data.savename="${prefix}${train_name}_on_Bloom2019_RMxYPS163_${suffix}" \
 data.test_data='${oc.env:DATA_DIR}/regression_data/new_latent/bloom2019_RMxYPS163_regression_std.feather' $extra_params
echo "${train_name} on Bloom2019_RMxYPS163"
