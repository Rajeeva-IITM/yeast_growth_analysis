#!/bin/bash

python src/compare_performance.py \
'run_path=${oc.env:RUN_DIR}/full/' \
'out_path=${oc.env:ROOT_DIR}/boosting/performance/Final/Full' \
'data_path=${oc.env:DATA_DIR}/full/new_latent' \
'data_paths={Bloom2013: ${data_path}/bloom2013_clf_3_pubchem.feather, Bloom2015: ${data_path}/bloom2015_clf_3_pubchem.feather}' \
'model_load_keys.prefix=Full_' 'model_load_keys.model_names=[Bloom2013, Bloom2015]' \
~data_paths.Bloom2019 ~data_paths.Bloom2019_BYxM22 ~data_paths.Bloom2019_RMxYPS163
