#!/usr/bin/env bash

# SLURM
#SBATCH --job-name=03_effect_of_priors
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --qos=short
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25
#SBATCH --mem-per-cpu=2G

# ENVIRONMENT #
# module load build-env/2020
module load anaconda3/2019.03
source $EBROOTANACONDA3/etc/profile.d/conda.sh
conda activate faps

srun python 03_analysis/05_simulations/01_simulate_fixed_family_sizes.py
