#!/usr/bin/env python3

"""
Tom Ellis, 1st October 2021

Script to infer mating events for each iteration of an MCMC output. 
"""
import numpy as np
import os

from amajusmating import mating

# FAPS objects and distance matrices are generated in a separate script.
exec(open('03_analysis/01_data_formatting/setup_FAPS_GPS.py').read())

# Input and output data
input_dir = "03_analysis/03_mcmc/01_mcmc_restrict_kurtosis/output/"
output_dir = os.path.dirname(os.path.abspath(__file__)) + '/output/'
os.makedirs(output_dir, exist_ok=True)

# Parameters for posterior simulations of mating
np.random.seed(87)
burnin = 1500

print("\nInferring mating events from each iteration of the MCMC chain, discarding the first {} iterations as burn-in.\n\n".format(burnin))
mating.mating_over_chains(
    data = am_data,
    input_dir = input_dir,
    output_dir = output_dir,
    burnin = burnin
    )

