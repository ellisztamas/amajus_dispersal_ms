# from warnings import warn
# from faps.alogsumexp import alogsumexp
from operator import pos
import numpy as np
import pandas as pd
import faps as fp
from glob import glob
from tqdm import tqdm
import os

def import_mcmc(folder, burnin):
    """
    Import files with MCMC output for A. majus mating parameters

    Glob a set of files with the output of amajusmating.mating.run_MCMC(), import
    each and remove the first rows of each as burn-in, then concatenate into a 
    single dataframe.

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
    # For this script we are only interested in the generalised-Gaussian part
    # of the dispersal kernel, so set `mixture` to 1 for all rows.
    posterior['mixture'] = 1

    return posterior

def simulate_mating(data, model, ndraws = 1000, max_distance = np.inf):
    """
    Posterior samples of mating events.

    Draw a sample of mating events from the posterior distributions of sibship
    stuctures and paternity, given observed genotype and dispersal data. Also 
    draw a sample of mating events based on dispersal only, to simulate what 
    would be 'expected' if mating was random, and based only on distance.
    Finally, merge the two datasets with GPS and phenotype data.

    Parameters
    ==========
    data: `faps_data` class object
        Data about the population.
    model: dict
        Dictionary of starting model parameters. Keys should be a subset of
        ['missing', 'shape', 'scale', 'mixture', 'assortment'] and values 
        floats giving initial values for those parameters, within appropriate
        boundaries.
    ndraws: int, optional
        Number of Monte Carlo draws to perform for sibship clustering. Defaults
        to 1000.
    max_distance: float, int
        Maximum distance from the mother a candidate may be. Candidates further than
        this value will have their posterior probability of paternity set to zero.
        This is equivalent to setting a threshold prior on distance.
    
    Returns
    =======
    A dictionary of two dataframes showing the output of faps.posterior_mating()
    for observed data ('obs'; based on genotypic and covariates) and expected
    data ('exp'; based on covariates only). This includes columns for flower
    colours and GPS.
    """
    # Update data with parameter values
    data.update_covariate_probs(model = model, max_distance = max_distance)
    # We need to set dispersal probabilities for mothers to zero
    # Otherwise the mother would be drawn many times when drawing from the dispersal
    # kernel only, and cause many invalid partitions.
    jx = [np.where(x == np.array(data.candidates))[0][0] for x in data.mothers]
    for i,j in zip(range(len(data.paternity)), jx):
        data.covariates['dispersal'][i, j] = -np.inf

    # Cluster into sibships, if not already done.
    data.sibship_clustering(ndraws = 1000, use_covariates = True)
    # Draw a posterior data set for observed data, and expected under dispersal only
    obs = fp.posterior_mating(data.sibships, covariates_only=False, ndraws=ndraws)
    exp = fp.posterior_mating(data.sibships, covariates_only=True, ndraws=ndraws)

    # Add data on distance between mother and father and flower colours.
    # Turn distance matrix into a data frame so we can use .loc on it.
    distance_df = pd.DataFrame({
        'mother'   : np.repeat(list(data.mothers), data.n_candidates),
        'father'   : np.tile(data.candidates, len(data.mothers)),
        'distance' : data.distances.flatten()
    })
    # Merge siring events with distances and phenotypes.
    # Observed mating events
    obs = obs.\
        merge(data.gps['Easting'], how='left', left_on='mother', right_index=True) .\
        merge(distance_df, how="left", on=['mother', 'father']).\
        merge(data.flower_colours, how="left", left_on="mother", right_index=True).\
        merge(data.flower_colours, how="left", left_on="father", right_index=True, suffixes = ['_mother', "_father"])
    # Random mating
    exp = exp.\
        merge(data.gps['Easting'], how='left', left_on='mother', right_index=True).\
        merge(distance_df, how="left", on=['mother', 'father']).\
        merge(data.flower_colours, how="left", left_on="mother", right_index=True).\
        merge(data.flower_colours, how="left", left_on="father", right_index=True, suffixes = ['_mother', "_father"])
        
    return {
        'obs' : obs,
        'exp' : exp
    }

def summarise_families(obs):
    """
    Summarise mating events and dispersal

    Parameters
    ==========
    obs: pandas.DataFrame
        Observed mating events given dispersal and genetic data.
        This should be the output from `amajus.mating.simulate_mating()`.
    
    Returns
    =======
    A list giving:
    
    * Number of mating events, excluding those with missing fathers.
    * Number of mating events with unsampled candidates. Note that if several 
        sires for a half-sib array are missing FAPS lumps these together, so the
        actual number could be higher.
    * Estimated number of 'orphans'; offspring with an unsampled father.
    * Number of half-sibships for which all fathers were unsampled.
    * Median pollen dispersal distance
    * Number of dispersal events > 100m
    * Number of dispersal events > 500m
    * Number of dispersal events > 1000m
    """
    return [
        np.sum(obs['father'] != "missing"),
        np.sum(obs['father'] == "missing"),
        obs.loc[obs['father'] == "missing"]['offspring_sired'].sum(),
        (obs.loc[obs['father'] == "missing"]['frequency'] == 1).sum(),
        np.median(obs['distance'].dropna()),
        np.sum(obs['distance'] > 100),
        np.sum(obs['distance'] > 500),
        np.sum(obs['distance'] > 1000),
    ]
    
def relative_fitness(obs, exp, column, geno, boundaries = None):
    """
    Relative fitness of genotypes across the whole hybrid zone, and in
    individual spatial bins.

    Relative fitness for the whole population is the absolute number of mating
    events sired by males of each genotype, divided by the maximum of those
    counts. This is run for the whole population. In addition, if a list of 
    Easting positions is given, the mothers are divided into spatial bins, and
    relative fitnesses of sires is calculated separately for each bin. This is
    done separately for observed and expected datasets

    Parameters
    ==========
    obs: pandas.DataFrame
        Observed mating events given dispersal and genetic data.
        This should be the output from `amajus.mating.simulate_mating()`.
    exp: pandas.DataFrame
        Expected mating events under dispersal only.
        This should be the output from `amajus.mating.simulate_mating()`.
    column: str
        Column name shared by data.mating['obs'] and data.mating['exp'] giving
        genotypes for which relative fitnesses should be estimated
    geno: list
        List of possible genotypes to identify in `column`.
    boundaries: list or int, optional
        List giving Easting positions to divide mothers into spatial bins. This
        is passed to the `bins` argument of `pandas.cut()`; see that function
        for details.
    
    Returns
    =======
    A dataframe giving relative fitnesses for each genotype overall and in each 
    spatial bin in observed and expected datasets. Columns indicate:
    
    - 'dataset': Observed vs expected datasets
    - 'bin': Label for the spatial bin. 'all' indicates rows for the whole hybrid zone
    - 'start': Western-most boundary of the bin
    - 'stop' Eastern-most boundary of the bin
    - 'n_sires': Number of sires
    - Subsequent columns : relative fitnesses of each genotype
    """
    # Fathers with missing data mess up relative fitnesses
    # Pull out only those rows with phenotype data
    obs = obs.query(column + ' in ' + str(geno))  
    exp = exp.query(column + ' in ' + str(geno))  
    
    # Divide mothers into spatial bins
    if boundaries:
        # Create a Series with an integer label for each bin. I used integers
        # instead of the standard output of pd.cut because it was easier to loop over.
        region_obs = pd.cut(obs['Easting'], bins = boundaries, labels = np.arange(1, len(boundaries)))
        region_exp = pd.cut(exp['Easting'], bins = boundaries, labels = np.arange(1, len(boundaries)))

    # Empty list to store relative fitnesses
    rel_fitness = []

    # Absolute fitness each genotype across the whole hybrid zone.
    # Sum mating events from sires of each genotype
    abs_obs = [(obs[column] == g).sum() for g in geno]
    abs_exp = [(exp[column] == g).sum() for g in geno]
    # Relative fitness of each genotype
    rel_obs = abs_obs / np.max(abs_obs)
    rel_exp = abs_exp / np.max(abs_exp)
    # Send to rel_fitness. NaNs are for bin start and stop, which aren't relevant here
    rel_fitness = rel_fitness + [['obs', 'all', np.nan, np.nan, np.sum(abs_obs)] + list(rel_obs)]
    rel_fitness = rel_fitness + [['exp', 'all', np.nan, np.nan, np.sum(abs_exp)] + list(rel_exp)]
    
    if boundaries:
        # Relative fitnesses in each spatial bin.
        for r in range(1, len(boundaries)):      
            # Absolute fitness each genotype (= sum mating events from sires of each genotype)
            abs_obs = [(obs.loc[region_obs == r][column] == g).sum() for g in geno]
            abs_exp = [(exp.loc[region_exp == r][column] == g).sum() for g in geno]
            # Relative fitness of each genotype
            rel_obs = abs_obs / np.max(abs_obs)
            rel_exp = abs_exp / np.max(abs_exp)
            # Send to rel_fitness, along with bin label, and start and stop for each bin
            rel_fitness = rel_fitness + [['obs', r, boundaries[r-1], boundaries[r], np.sum(abs_obs)] + list(rel_obs)]
            rel_fitness = rel_fitness + [['exp', r, boundaries[r-1], boundaries[r], np.sum(abs_exp)] + list(rel_exp)]
    
    return pd.DataFrame(rel_fitness, columns=['dataset', 'bin', 'start', 'stop', 'n_sires'] + geno)

def assortative_mating(obs, exp, col_x, col_y, boundaries = None):
    """
    Assortative mating between genotypes across the hybrid zone, and within
    spatial bins.

    Assortment probabilities for the whole population are calculated as the mean
    number of mating events between individuals of the same genotype. This is
    run for the whole population. In addition, if a list of Easting positions is
     given, the mothers are divided into spatial bins, and assortment 
    probabilities are calculated separately for each bin. This is done
    separately for observed and expected datasets.
    
    Parameters
    ==========
    obs: pandas.DataFrame
        Observed mating events given dispersal and genetic data.
        This should be the output from `amajus.mating.simulate_mating()`.
    exp: pandas.DataFrame
        Expected mating events under dispersal only.
        This should be the output from `amajus.mating.simulate_mating()`.
    col_x, col_y: str
        The pair of column names shared by data.mating['obs'] and
        data.mating['exp'] giving genotypes for which assortment probabilities
        should be estimated.
    boundaries: list or int, optional
        List giving Easting positions to divide mothers into spatial bins. This
        is passed to the `bins` argument of `pandas.cut()`; see that function
        for details.

    Returns
    =======
    A dataframe giving relative fitnesses for each genotype overall and in each 
    spatial bin in observed and expected datasets. Columns indicate:
    
    - 'dataset': Observed vs expected datasets
    - 'bin': Label for the spatial bin. 'all' indicates rows for the whole
        hybrid zone
    - 'start': Western-most boundary of the bin
    - 'stop' Eastern-most boundary of the bin
    - 'n_sires': Number of sires
    - 'assortment : Mean mating events between individuals of the same genotype.
    
    """
    # Remove rows with missing data for one or both parents
    obs = obs.copy().loc[obs[col_x].notna() & obs[col_y].notna()]
    exp = exp.copy().loc[exp[col_x].notna() & exp[col_y].notna()]
    # Create a column stating whether genotypes match or not
    obs['match'] = obs[col_x] == obs[col_y]
    exp['match'] = exp[col_x] == exp[col_y]

    # Divide mothers into spatial bins
    if boundaries:
        # Create a Series with an integer label for each bin. I used integers
        # instead of the standard output of pd.cut because it was easier to loop
        # over.
        region_obs = pd.cut(obs['Easting'], bins = boundaries, labels = np.arange(1, len(boundaries)))
        region_exp = pd.cut(exp['Easting'], bins = boundaries, labels = np.arange(1, len(boundaries)))
    
    # Empty list to store assortment probabilities
    assortment = []

    # Assortative mating for the whole hybrid zone
    ass_obs = obs['match'].mean()
    ass_exp = exp['match'].mean()
    # Send to assortment. NaNs are for bin start and stop, which aren't relevant here
    assortment = assortment + [['obs', 'all', np.nan, np.nan, obs.shape[0], ass_obs]]
    assortment = assortment + [['exp', 'all', np.nan, np.nan, exp.shape[0], ass_exp]]

    # Assortative mating in each spatial bin.
    if boundaries:
        for r in range(1, len(boundaries)):      
            # Assortment probabilities in each dataset
            ass_obs = obs.loc[region_obs == r]['match'].mean()
            ass_exp = exp.loc[region_exp == r]['match'].mean()
            # Send to assortment, along with bin label, start and stop, and sample sizes for each bin
            assortment = assortment + [['obs', r, boundaries[r-1], boundaries[r], (region_obs == r).sum(), ass_obs]]
            assortment = assortment + [['exp', r, boundaries[r-1], boundaries[r], (region_exp == r).sum(), ass_exp]]

    return pd.DataFrame(
        assortment,
        columns=['dataset', 'bin', 'start', 'stop', 'n_sires', 'assortment']
        )

def mating_over_chains(data, folder, boundaries, burnin = 500, ndraws = 1000):
    """
    Summarise number of mating events, relative fitness and assortment for each
    iteration of a set of MCMC files.

    Parameters
    ==========
    data: `faps_data` class object
        Data about the population.
    folder: str
        Path to a directory containing one or more files from 
        amajusmating.mating.run_MCMC ending with the suffix `.out`
    boundaries: list or int, optional
        List giving Easting positions to divide mothers into spatial bins. This
        is passed to the `bins` argument of `pandas.cut()`; see that function
        for details.
    burnin: int
        Integer number of rows to remove from each chain as burnin
    ndraws: int, optional
        Number of Monte Carlo draws to perform for sibship clustering. Defaults
        to 1000.

    Returns
    =======
    Saves CSV files for the outputs of
    
    * mating.summarise_mating()
    * mating.relative_fitness() for phenotype data (full red, hybrid, yellow),
        plus Rosea and Sulfurea genotypes
    * assortment() for phenotypes (full red, hybrid, yellow)

    Also saves all mating events as `sires.csv`
    """
    # Import MCMC results 
    mcmc = import_mcmc(folder, burnin = burnin)
    # Empty dictionaries to store the output of each iteration
    sires      = {}
    summarise  = {}
    phenotype  = {}
    rosea      = {}
    sulf       = {}
    assortment = {}

    # Loop over steps in the MCMC chain and get siring events for each.
    for i in tqdm(mcmc.index):
        
        model = mcmc.loc[i]
        # Simulate mating events, and include GPS and phentoype information about mother and sire
        data.mating = simulate_mating(data, model, ndraws=ndraws)

        # Full list of mothers, fathers and distances between them.
        sires[i] = data.mating['obs'][['mother','father','distance']]
        # Summarise mating events and missing fathers.
        summarise[i] = summarise_families(data.mating['obs'])
        # Relative fitness of full red, hybrid and yellow plants
        phenotype[i] = relative_fitness(
            obs = data.mating['obs'],
            exp = data.mating['exp'],
            column = "simple_colour_father",
            geno = ["FR", "hybrid", "Ye"],
            boundaries  = boundaries
        )
        # Relative fitness of Rosea genotypes
        rosea[i] = relative_fitness(
            obs = data.mating['obs'],
            exp = data.mating['exp'],
            column = "rosea_father",
            geno = ["R/R", "R/r", "r/r"],
            boundaries  = boundaries
        )
        # Relative fitness of Sulfurea genotypes
        sulf[i] = relative_fitness(
            obs = data.mating['obs'],
            exp = data.mating['exp'],
            column = "sulfurea_father",
            geno = ["S/+", "s/s"],
            boundaries  = boundaries
        )
        # Assortative mating.
        assortment[i] = assortative_mating(
            obs = data.mating['obs'],
            exp = data.mating['exp'],
            col_x = 'simple_colour_mother',
            col_y = 'simple_colour_father',
            boundaries  = boundaries
        )

    # Concatenate output for each iteration into single dataframes
    sires      = pd.concat(sires)
    phenotype  = pd.concat(phenotype)
    rosea      = pd.concat(rosea)
    sulf       = pd.concat(sulf)
    assortment = pd.concat(assortment)
    # Concatenate `summarise`, and also add column names.
    summarise = pd.DataFrame(summarise).T
    summarise.columns = [
        'n_mating_events',
        'missing_dads',
        'orphans',
        'empty_sibships',
        'median_dispersal',
        'dispersal>100m',
        'dispersal>500m',
        'dispersal>1000m'
    ]

    sires.to_csv(folder + "sires.csv")
    summarise.to_csv(folder + "summarise_mating.csv")
    phenotype.to_csv(folder + "phenotype_relative_fitness.csv")
    rosea.to_csv(folder + "rosea_relative_fitness.csv")
    sulf.to_csv(folder + "sulf_relative_fitness.csv")
    assortment.to_csv(folder + "assortment.csv")
