#!/bin/bash

#SBATCH --job-name=knit_rmd           # Job name
#SBATCH --output=logs/slurm/knit_rmd_%j.out # Output file
#SBATCH --error=logs/slurm/knit_rmd_%j.err  # Error file
#SBATCH --time=01:00:00               # Time limit (1 hour)
#SBATCH --mem=32G                     # Memory limit (32 GB)
#SBATCH --cpus-per-task=1             # Number of CPUs

module load R # Load R module if needed

RMD_FILE_PATH=$1 # Get the first argument passed to the script

# Run the R script and redirect its output and error to specific files
Rscript knit_rmd.R $RMD_FILE_PATH > logs/r/r_knit_rmd_output.log 2> logs/r/r_knit_rmd_error.log