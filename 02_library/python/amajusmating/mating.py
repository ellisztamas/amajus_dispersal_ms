import numpy as np
import pandas as pd
from glob import glob
from tqdm import tqdm
import os
from time import strftime
import platform
import weighted

import faps as fp

def import_mcmc(folder, burnin, mixture_is_1 = True):
    """
    Import and trim multiple files with MCMC output for A. majus mating
    parameters.

    This looks in a folder with one or more files containing the output of 
    amajusmating.mating.run_MCMC(), which end with `.out`, usually from separate
    chains from the same analysis. It imports these, and removes entries from 
    the first entries as burnin, then concatenates remaining steps into a single
    dataframe.

    Parameters
    ==========
    folder: str
        Path to a directory containing one or more files from 
        amajusmating.mating.run_MCMC ending with the suffix `.out`
    burnin: int
        Integer number of rows to remove from each chain as burnin

    Returns
    =======
    A single dataframe concatenating data from the imported chains. Since we are
    interested in dispersal parameters, the `mixture` parameter is set to 1 for
    all rows.
    """
    # Import MCMC chains
    chains = glob(folder+"/*out")

    posterior = {}
    for chain in chains:
        k = os.path.basename(chain)
        posterior[k] = pd.read_csv(chain, sep="\t").loc[lambda x: x.iter >= burnin]
    posterior = pd.concat(posterior).reset_index()
    if mixture_is_1:
        # If we are only interested in the generalised-Gaussian part
        # of the dispersal kernel, so set `mixture` to 1 for all rows.
        posterior['mixture'] = 1

    return posterior

def get_mating_events(data, model, max_distance = np.inf):
    """
    Generate a dataframe of mating events between mothers and plausible 
    candidate fathers, along with posterior support for each event having
    happened.

    This calls `faps.summarise_sires()` using observed genotype and dispersal
    data. For each maternal family, this pulls out a list of the plausible
    fathers associated with each possible partition structure, and the posterior
    probability of that structure. The probability of a mating event between the
    mother and an individual candidate father is then the sum of those 
    probabilities for each partition structure in which he appears.

    Likewise the number of offspring sired is the weighted mean number of
    offspring assigned to that father across partition structures, weighted by 
    the posterior probabilty of each partition.

    Parameters
    ==========
    data: `faps_data` class object
        Data about the population.
    model: dict
        Dictionary of starting model parameters. Keys should be a subset of
        ['missing', 'shape', 'scale', 'mixture', 'assortment'] and values 
        floats giving initial values for those parameters, within appropriate
        boundaries.
    max_distance: float, int
        Maximum distance from the mother a candidate may be. Candidates further
        than this value will have their posterior probability of paternity set 
        to zero. This is equivalent to setting a threshold prior on distance.
        Not used.
    
    Returns
    =======
    A dataframe of plausible mating events, showing the mother, father, 
    posterior support for that mating event, posterior mean number of offspring
    generated, as well as GPS and flower-colour information
    """
    # Update data with parameter values
    data.update_covariate_probs(model = model, max_distance = max_distance)
    # We need to set dispersal probabilities for mothers to zero
    # Otherwise the mother would be drawn many times when drawing from the dispersal
    # kernel only, and cause many invalid partitions.
    jx = [np.where(x == np.array(data.candidates))[0][0] for x in data.mothers]
    for i,j in zip(range(len(data.paternity)), jx):
        data.covariates['dispersal'][i, j] = -np.inf
    # Turn distance matrix into a data frame so we can use .loc on it.
    distance_df = pd.DataFrame({
        'mother'   : np.repeat(list(data.mothers), data.n_candidates),
        'father'   : np.tile(data.candidates, len(data.mothers)),
        'distance' : data.distances.flatten()
    })
    # Cluster into sibships, if not already done.
    data.sibship_clustering(ndraws = 1000, use_covariates = True)

    # Siring events based on genetics and dispersal    
    mating_events = fp.summarise_sires(data.sibships)
    # Merge siring events with distances and phenotypes.
    mating_events = mating_events.\
        merge(data.gps['Easting'], how='left', left_on='mother', right_index=True) .\
        merge(data.gps['Easting'], how='left', left_on='father', right_index=True, suffixes = ['_mother', "_father"]).\
        merge(distance_df, how="left", on=['mother', 'father'])
    
    return mating_events

def summarise_families(mating_events):
    """
    A function to create a summary of information about mating patterns in each 
    iteration of an MCMC chain. See `Returns` below for the data described.

    Parameters
    ==========
    mating_events: pandas.DataFrame
        Observed mating events given dispersal and genetic data.
        This should be the output from `amajusmating.simulate_mating()`.
    
    Returns
    =======
    A list giving:
    
    * Number of mating events, excluding those with missing fathers.
    * Number of mating events with unsampled candidates. Note that if several 
        sires for a half-sib array are missing FAPS lumps these together, so the
        actual number could be higher.
    * Estimated number of 'orphans'; offspring with an unsampled father.
    * Median pollen dispersal distance
    * Number of dispersal events > 100m
    * Number of dispersal events > 500m
    * Number of dispersal events > 1000m
    """
    return [
        mating_events['prob'].loc[~mating_events['father'].isna()].sum(), # Number of mating events with sampled fathers
        mating_events['prob'].loc[mating_events['father'].isna()].sum(), # Mating events with unsampled fathers
        mating_events['offspring'].loc[mating_events['father'].isna()].sum(), # Number of orphans
        weighted.median(
            data    = mating_events['distance'].loc[~mating_events['father'].isna()],
            weights = mating_events['prob'].loc[~mating_events['father'].isna()]
            ),
        mating_events['prob'].loc[( ~mating_events['father'].isna() ) & (mating_events['distance'] > 100)].sum(),
        mating_events['prob'].loc[( ~mating_events['father'].isna() ) & (mating_events['distance'] > 500)].sum(),
        mating_events['prob'].loc[( ~mating_events['father'].isna() ) & (mating_events['distance'] > 1000)].sum(),
    ]

def mating_over_chains(data, input_dir, output_dir, burnin = 500, ndraws = 1000):
    """
    For a folder of MCMC results, this imports the data, infers mating events
    for each iteration, summarises the results and saves the output.

    Parameters
    ==========
    data: `faps_data` class object
        Data about the population.
    input_dir: str
        Path to a directory containing one or more files from 
        amajusmating.mating.run_MCMC ending with the suffix `.out`
    output_dir: str
        Path to the directory where output should be saved.
    burnin: int
        Integer number of rows to remove from each chain as burnin
    ndraws: int, optional
        Number of Monte Carlo draws to perform for sibship clustering. Defaults
        to 1000.

    Returns
    =======
    Saves CSV files for the outputs of `get_mating_events()` and 
    `summarise_mating()`. Also saves a log file.
    """
    # Set up log file
    log_file = open(output_dir + "log", 'w')
    log_file.write(
        'Analysis of relative fitness and assortative mating in A. majus begun {} using FAPS {} in Python {}.\n\n'.format(
            strftime("%Y-%m-%d %H:%M:%S"),
            fp.__version__,
            platform.python_version()
            )
        )
    log_file.write('Using MCMC results from:\n{}\n\n'.format(input_dir))
    log_file.write('Saving output to:\n{}\n\n'.format(output_dir))
    log_file.write('Discarding the first {} iterations from each chain as burnin.'.format(burnin))
    log_file.close()

    # Import MCMC results 
    mcmc = import_mcmc(input_dir, burnin = burnin)
    # Empty dictionaries to store the output of each iteration
    sires      = {}
    summarise  = {}

    # Loop over steps in the MCMC chain and get siring events for each.
    for i in tqdm(mcmc.index):
        
        model = mcmc.loc[i]
        # Simulate mating events, and include GPS and phentoype information about mother and sire
        data.mating = get_mating_events(data, model)
        data.mating.insert(loc=0, column='iter', value=i)

        # Full list of mothers, fathers and the probabilities for each mating event.
        sires[i] = data.mating
        # Summarise mating events and missing fathers.
        summarise[i] = summarise_families(data.mating)
        # Relative fitness of full red, hybrid and yellow plants

    # Concatenate output for each iteration into single dataframes
    sires      = pd.concat(sires)
    # Concatenate `summarise`, and also add column names.
    summarise = pd.DataFrame(summarise).T
    summarise.columns = [
        'n_mating_events',
        'missing_dads',
        'orphans',
        'median_dispersal',
        'dispersal>100m',
        'dispersal>500m',
        'dispersal>1000m'
    ]

    sires.to_csv(    output_dir + "mating_events_over_chains.csv", float_format='%.3f', index=False)
    summarise.to_csv(output_dir + "summarise_mating.csv", float_format='%.3f', index=False)