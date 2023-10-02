#!/usr/bin/env bash

# SLURM
#SBATCH --mem=10GB
#SBATCH --job-name=amajus_mating_events
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --qos=medium
#SBATCH --time=1-00:00:00

# ENVIRONMENT #
# module load build-env/2020
module load anaconda3/2019.03
source $EBROOTANACONDA3/etc/profile.d/conda.sh
conda activate faps

srun python 03_analysis/05_simulations/01_simulate_mating_events.py
