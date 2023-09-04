#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=2g
#SBATCH -t 5-00:00:00
#SBATCH --chdir=/proj/raulab/users/brian/r_projects/cc_composition  
#SBATCH --output=./slurm/logs/


# Set up job environment:
module purge   # clear any inherited modules
module load r/4.2.1

# Run Rscript in a clean R instance, output a logfile
Rscript --vanilla --verbose ./scripts/utilities/r_libs.R > ./slurm/logs/slurm-${SLURM_JOBID}.Rout 2>&1 

# append logfile to this scripts logfile
cat ./slurm/logs/slurm-${SLURM_JOBID}.Rout >> ./slurm/logs/slurm-${SLURM_JOBID}.out

# remove Rout log
rm ./slurm/logs/slurm-${SLURM_JOBID}.Rout
