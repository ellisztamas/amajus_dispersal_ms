#!/usr/bin/env bash

# SLURM
#SBATCH --job-name=mcmc_no_mixture
#SBATCH --mem=10GB
#SBATCH --output=slurm/mcmc_no_mixture.%J.out
#SBATCH --error=slurm/mcmc_no_mixture.%J.err
#SBATCH --qos=long
#SBATCH --time=3-00:00:00
#SBATCH --array=1-4

# ENVIRONMENT #
module load build-env/2020
module load r/3.5.1-foss-2018b
module load anaconda3/2019.03
source $EBROOTANACONDA3/etc/profile.d/conda.sh
conda activate faps

srun python 03_analysis/03_mcmc/04_no_mixture/run_MCMC.py -i $SLURM_ARRAY_TASK_ID
