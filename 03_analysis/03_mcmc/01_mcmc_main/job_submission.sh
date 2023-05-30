#!/usr/bin/env bash

# SLURM
#SBATCH --job-name=mcmc_main
#SBATCH --mem=20GB
#SBATCH --output=./03_analysis/03_mcmc/mcmc_main.%J.out
#SBATCH --error=./03_analysis/03_mcmc/mcmc_main.%J.err
#SBATCH --qos=medium
#SBATCH --time=3-00:00:00
#SBATCH --array=0-3

# ENVIRONMENT #
module load build-env/2020
module load r/3.5.1-foss-2018b
module load anaconda3/2019.03
source $EBROOTANACONDA3/etc/profile.d/conda.sh
conda activate faps

srun python /03_analysis/03_mcmc/_01_mcmc_main/run_MCMC.py -i $SLURM_ARRAY_TASK_ID
