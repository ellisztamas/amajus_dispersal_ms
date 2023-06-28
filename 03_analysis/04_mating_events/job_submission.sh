#!/usr/bin/env bash

# SLURM
#SBATCH --job-name=mating_events
#SBATCH --mem=10GB
#SBATCH --qos=medium
#SBATCH --output=./03_analysis/04_mating_events/log.txt
#SBATCH --time=1-00:00:00
#SBATCH --array=0,3

# ENVIRONMENT #
module load build-env/2020
module load anaconda3/2019.03
source $EBROOTANACONDA3/etc/profile.d/conda.sh
conda activate faps

FILES=(./03_analysis/03_mcmc/**)

srun python ./03_analysis/04_mating_events/get_mating_events.py -i ${FILES[$SLURM_ARRAY_TASK_ID]}
