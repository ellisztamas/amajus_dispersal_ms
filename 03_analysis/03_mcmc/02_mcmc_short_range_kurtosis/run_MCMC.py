#!/usr/bin/env python3

"""
Tom Ellis, 26th January 2021

Script to run joint analysis of paternity, sibships and dispersal using
priors that allow kurtosis but at short scales
"""

import numpy as np
import os
from scipy.stats import beta
from scipy.stats import gamma

from amajusmating import mcmc

# FAPS objects and distance matrices are generated in a separate script.
exec(open('03_analysis/01_data_formatting/setup_FAPS_GPS.py').read())

# INITIALISE THE MODEL
nreps = 4000 # Total number of iterations to run
thin  = 10 # How often to write samples.
max_distance = np.inf # set a maximum dispersal distance
# output_dir = "005.results/004_mcmc_restrict_kurtosis/output/"
output_dir = os.path.dirname(os.path.abspath(__file__))+'/output/'
os.makedirs(output_dir, exist_ok=True)

np.random.seed(1246)

# Parameters for posterior simulations of mating
ndraws = 1000
burnin = 500
spatial_bins = [-827, -70, 70, 244, 270, 560]

# PRIORS
priors = (lambda x : {
    'missing' : beta.pdf(x['missing'], a=3,   b=15),
    'mixture' : beta.pdf(x['mixture'], a=1.1, b=1.1),
    'shape'   : gamma.pdf(x['shape'],   a=2,  scale = 1/2),
    'scale'   : gamma.pdf(x['scale'],   a=2,   scale = 3)
})

# Proposed values are a Gaussian peturbation away from the previous values.
# This is controlled by the sigma of the gaussian, which is defined for each variable
proposal_sigma = {
    'missing' : 0.0,
    'shape'  : 0.05,
    'scale'  : 2,
    'mixture' : 0.025,
}

for i in [1,2,3,4]:
    mcmc.run_MCMC(
        data= am_data,
        initial_parameters = {
            'missing' : [0.32,0.32,0.32,0.32] [i-1],
            'shape'   : [   2, 0.5, 0.2,   1] [i-1],
            'scale'   : [  70,  40, 100,  10] [i-1],
            'mixture' : [0.99, 0.4, 0.8, 0.6] [i-1]
        },
        proposal_sigma = proposal_sigma,
        priors = priors,
        thin=thin,
        nreps=nreps,
        output_dir = output_dir,
        chain_name = 'chain' + str(i),
        max_distance = max_distance
        )