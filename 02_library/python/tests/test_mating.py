"""
format_mcmc
    do this in R so you can also check convergence etc
    give it a folder name, burnin, output sample size
    trim burnin and thin
    concatenate
observed_sires 
    draw dads using genetics and dispersal for a given model
    return a df
random_sires
    to draw dads based on dispersal only, for a given model
    add max_distance
    add update covariates function, and add model into here.
    return a df
posterior_mating
    loop over mcmc output
    convert mcmc line into a model as a dict
    check that those parameters are valid.
    call observed_sires and random_sires on each
    return df for observed and random sires.
"""
import numpy as np
import pandas as pd

from amajusmating import mating


# FAPS objects and distance matrices are generated in a separate script.
exec(open('03_analysis/01_data_formatting/setup_FAPS_GPS.py').read())

model = {
    'missing' : 0.29,
    "mixture" : 0.8,
    "shape"   : 0.6,
    "scale"   : 25
}

def test_observed_sires():
    obs = mating.observed_sires(data= am_data, model = model)
    assert isinstance(obs, pd.core.frame.DataFrame)