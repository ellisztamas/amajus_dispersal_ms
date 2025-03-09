#!/usr/bin/env bash

# SLURM
#SBATCH --job-name=mcmc_no_mixture
#SBATCH --mem=10GB
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --qos=short
#SBATCH --time=8:00:00
#SBATCH --array=1-4

# ENVIRONMENT #
source activate amajus_mating_ms

srun python 03_analysis/03_mcmc/02_no_mixture/run_MCMC.py -i $SLURM_ARRAY_TASK_ID
