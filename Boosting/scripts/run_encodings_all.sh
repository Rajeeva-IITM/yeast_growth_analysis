#!/bin/bash

# All three gene encodings, in increasing order of cost. Run from inside a `pixi shell`,
# from the project root. Defaults to sigma_0.5; pass 1.0 for the other threshold family.
#
#     ./scripts/run_encodings_all.sh
#     ./scripts/run_encodings_all.sh 1.0
#
# max and sum are 6014 gene columns each; count is 18042 and will take far longer than the
# other two combined. Run them in this order so the baseline numbers land first.

alpha="${1:-0.5}"

./scripts/run_encodings_max.sh "${alpha}"
./scripts/run_encodings_sum.sh "${alpha}"
./scripts/run_encodings_count.sh "${alpha}"

echo "All encodings Done"
