#!/bin/bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -t 48:00:00
#SBATCH --mem=1g
#SBATCH -n 1

module load r r/4.2.1

# Run the R script with the model index and chunk index as arguments
Rscript scripts/benchmarking/deseq/parallel_models/09182023/single_model.R $1 $2
