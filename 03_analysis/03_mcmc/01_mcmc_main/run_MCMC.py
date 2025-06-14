#!/usr/bin/env python3

"""
Script to run joint analysis of paternity, sibships and dispersal by Metropolis-
Hastings MCMC.

Tom Ellis, 3rd April 2023
"""
import numpy as np
import os
from scipy.stats import beta
from scipy.stats import gamma
from scipy.stats import lognorm
import argparse

from amajusmating import mcmc

parser = argparse.ArgumentParser()
parser.add_argument('-i', required = True, type=int)
args = parser.parse_args()

# FAPS objects and distance matrices are generated in a separate script.
exec(open('03_analysis/01_data_formatting/setup_FAPS_GPS.py').read())

# INITIALISE THE MODEL
nreps = 25500 # Total number of iterations to run
thin  = 100 # How often to write samples.
max_distance = np.inf # set a maximum dispersal distance
output_dir = os.path.dirname(os.path.abspath(__file__))+'/output/'
os.makedirs(output_dir, exist_ok=True)

np.random.seed(87)

# PRIORS
priors = (lambda x : {
    'missing' : beta.pdf(x['missing'], a=3,   b=15),
    'mixture' : beta.pdf(x['mixture'], a=1.1, b=1.1),
    'shape'   : lognorm.pdf(x['shape'], scale = 0.5, s=0.5),
    'scale'   : gamma.pdf( x['scale'], a=2,   scale = 100)
})

# Proposed values are a Gaussian peturbation away from the previous values.
# This is controlled by the sigma of the gaussian, which is defined for each variable
proposal_sigma = {
    'missing' : 0.0,
    'shape'  : 0.05,
    'scale'  : 2,
    'mixture' : 0.025,
}

print("\nBeginning MCMC.\n\n")

i = args.i
mcmc.run_MCMC(
    data= am_data,
    initial_parameters = {
        'missing' : [ 0.5, 0.5, 0.5, 0.5] [i-1],
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
