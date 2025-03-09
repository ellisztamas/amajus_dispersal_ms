#!/usr/bin/env bash

# SLURM
#SBATCH --job-name=mating_events
#SBATCH --mem=10GB
#SBATCH --qos=short
#SBATCH --output=slurm/%x-%a.out
#SBATCH --error=slurm/%x-%a.err
#SBATCH --time=8:00:00
#SBATCH --array=0

# ENVIRONMENT #
source activate amajus_mating_ms

FILES=(03_analysis/03_mcmc/**)

srun python ./03_analysis/04_mating_events/get_mating_events.py -i ${FILES[$SLURM_ARRAY_TASK_ID]}
