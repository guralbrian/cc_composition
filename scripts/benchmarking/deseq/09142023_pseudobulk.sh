#!/bin/bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -t 48:00:00
#SBATCH --mem=8g
#SBATCH -n 1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bgural@email.unc.edu

module load r r/4.2.1

# Run the R script with the start and end index as arguments
Rscript scripts/benchmarking/deseq/09142023_pseudobulk.R 