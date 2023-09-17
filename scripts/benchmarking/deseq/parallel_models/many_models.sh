#!/bin/bash

#!/bin/bash
# Loop over models
for model in $(seq $1 $2); do
  # Loop over chunks
  for chunk in $(seq 1 20); do  # Here 20 should be ceil(10000 / 500)
    sbatch scripts/benchmarking/deseq/parallel_models/single_model.sh $model $chunk
  done
done
