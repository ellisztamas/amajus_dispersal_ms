#!/usr/bin/env bash

# SLURM
#SBATCH --job-name=mcmc_main
#SBATCH --mem=20GB
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --qos=long
#SBATCH --time=3-00:00:00
#SBATCH --array=1-4

# ENVIRONMENT #
module load build-env/f2022
module load anaconda3/2023.03
source ~/.bashrc
conda activate faps

srun python 03_analysis/03_mcmc/01_mcmc_main/run_MCMC.py -i $SLURM_ARRAY_TASK_ID
