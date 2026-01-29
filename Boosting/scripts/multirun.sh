#!/bin/bash

# cd "/home/rajeeva/Project/boosting/" || exit
# There's a better way to do it

for seed in {1..5}
do
    export SEED=${seed}
    source ./scripts/run_script.sh
done
