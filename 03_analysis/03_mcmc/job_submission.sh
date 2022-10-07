#!/usr/bin/env bash

# SLURM
#SBATCH --mem=40GB
#SBATCH --output=./03_analysis/03_mcmc/log.txt
#SBATCH --qos=medium
#SBATCH --time=24:00:00
#SBATCH --array=0-4

# ENVIRONMENT #
module load build-env/2020
module load r/3.5.1-foss-2018b
module load anaconda3/2019.03
source $EBROOTANACONDA3/etc/profile.d/conda.sh
conda activate faps

FILES=(./03_analysis/03_mcmc/**/run_MCMC.py)

srun python ${FILES[$SLURM_ARRAY_TASK_ID]}