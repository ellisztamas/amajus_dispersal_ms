#!/usr/bin/env bash

# SLURM
#SBATCH --mem=10GB
#SBATCH --qos=medium
#SBATCH --time=00:10:00
#SBATCH --array=0-5

# ENVIRONMENT #
# module load build-env/2020
module load anaconda3/2019.03
source $EBROOTANACONDA3/etc/profile.d/conda.sh
conda activate faps

FILES=(./03_analysis/03_mcmc/**/output/)

srun python ./03_analysis/04_mating_events/get_mating_events.py ${FILES[$SLURM_ARRAY_TASK_ID]}