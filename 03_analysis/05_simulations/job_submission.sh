#!/usr/bin/env bash

# SLURM
#SBATCH --mem=10GB
#SBATCH --job-name=amajus_simulations
#SBATCH --qos=medium
#SBATCH --output=slurm/amajus_simulations.out
#SBATCH --error=/amajus_simulations.err
#SBATCH --time=1-00:00:00

# ENVIRONMENT #
# module load build-env/2020
module load anaconda3/2019.03
source $EBROOTANACONDA3/etc/profile.d/conda.sh
conda activate faps

srun python 03_analysis/05_simulations/simulations.py
