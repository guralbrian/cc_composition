#!/bin/bash

# Loop over models 
# should be index positions for list of models within single_model.R
# i.e. there are 11 models, so run with "1 11" as trailing args
for model in $(seq $1 $2); do
  # Loop over chunks
  for chunk in $(seq 1 4); do  # Here 20 should be ceil(2000 / 500)
    sbatch scripts/benchmarking/deseq/parallel_models/09182023/single_model.sh $model $chunk
  done
done
