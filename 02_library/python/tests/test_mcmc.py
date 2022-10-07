import numpy as np
import os
import pandas as pd
from scipy.stats import beta
from scipy.stats import gamma as gma

from amajusmating import mcmc

# FAPS objects and distance matrices are generated in a separate script.
exec(open('03_analysis/01_data_formatting/setup_FAPS_GPS.py').read())

# INITIALISE THE MODEL
np.random.seed(1246)

# Dictionary listing starting values.
initial_parameters = {
    'missing' : 0.15, # proportion missing fathers
    'shape'  : 1,
    'scale'  : 10,
    'mixture' : 0.8
}

# Proposed values are a Gaussian peturbation away from the previous values.
# This is controlled by the sigma of the gaussian, which is defined for each variable
proposal_sigma = {
    'missing' : 0.025,
    'shape'  : 0.05,
    'scale'  : 2,
    'mixture' : 0.025,
}

# PRIORS
priors = (lambda x : {
    'missing' : beta.pdf(x['missing'], a=3,   b=15),
    'mixture' : beta.pdf(x['mixture'], a=1.1, b=1.1),
    'shape'   : gma.pdf(x['shape'],   a=10,  scale = 1/5),
    'scale'   : gma.pdf(x['scale'],   a=6,   scale = 50)
})

def test_mcmc():
    folder = os.path.dirname(os.path.abspath(__file__))
    file = "/mcmc_test_chain"

    chain = mcmc.run_MCMC(
        data= am_data,
        initial_parameters = initial_parameters,
        proposal_sigma = proposal_sigma,
        priors = priors,
        thin=1,
        nreps=3,
        output_dir = folder,
        chain_name = file,
        max_distance = np.inf
        )

    assert os.path.exists(folder + file + ".out")
    assert os.path.exists(folder + file + ".log")

    dat = pd.read_csv("002.library/python/tests/mcmc_test_chain.out", sep="\t")
    assert list(dat.keys()) == ['iter', 'hours', 'log_posterior', 'log_prior', 'loglik', 'missing', 'mixture', 'scale', 'shape']

    os.remove(folder + file + ".out")
    os.remove(folder + file + ".log")
