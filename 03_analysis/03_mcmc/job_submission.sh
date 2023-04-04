#!/usr/bin/env bash

# SLURM
#SBATCH --job-name=mcmc_jobs
#SBATCH --mem=10GB
#SBATCH --output=./03_analysis/03_mcmc/mcmc_jobs.%J.out
#SBATCH --error=./03_analysis/03_mcmc/mcmc_jobs.%J.err
#SBATCH --qos=long
#SBATCH --time=4-00:00:00
#SBATCH --array=0-4

# ENVIRONMENT #
module load build-env/2020
module load r/3.5.1-foss-2018b
module load anaconda3/2019.03
source $EBROOTANACONDA3/etc/profile.d/conda.sh
conda activate faps

FILES=(./03_analysis/03_mcmc/**/run_MCMC.py)

srun python ${FILES[$SLURM_ARRAY_TASK_ID]}