import numpy as np
import pandas as pd
import faps as fp
import warnings


def draw_sires(data, model, mating_events):
    """
    Draw vectors of candidates for each mother.
    
    For each maternal family, this draws integer indices of a sample of pollen
    donors to mate with, in proportion to dispersal probabilities given in 
    `model`.

    Parameters
    ==========
    data: `faps_data` class object
        Data about the population.
    model: dict
        Dictionary of starting model parameters. Keys should be a subset of
        ['missing', 'shape', 'scale', 'mixture', 'assortment'] and values 
        floats giving initial values for those parameters, within appropriate
        boundaries.
    mating_events: pandas.DataFrame
        Dataframe listing mating events. Usually this is the output of
        `amajusmating.mating.mating_over_chains()`. At a minimum it should list
        columns 'mother', 'father', 'offspring' (number of offspring in each
        family).
    
    Returns
    =======
    A list of integer indexes for each sire drawn, with one element for each
    row in `mating_events`.
    
    """
    # Include mating probabilities in faps_data object
    data.update_covariate_probs(model)
    
    # copy rows of dispersal probabilities to match the number of mating events
    ix = [ np.where(x == np.array(data.mothers))[0][0] for x in mating_events['mother'] ]

    assert len(ix) == mating_events.shape[0]

    # Zip together family sizes, and mating probabilities to iterate over later
    probs = np.exp(data.covariates['dispersal'][ix])
    
    # For each mother, draw a sample of pollen donors in proportion to their distance
    sires = []
    for p in probs:
        these_sires = np.random.choice(
            a = range(p.shape[0]), 
            size = 1,
            replace = False,
            p = p
        )
        sires = sires + [these_sires]
    sires = [x for y in sires for x in y]

    return [ data.candidates[i] for i in sires ]

def simulate_genotypes(mating_events, genotypes, mu):
    """
    This simulates families of full siblings of the same sizes as in 
    the real data set.
    
    Parameters
    ==========
    mating_events: pandas.DataFrame
        Dataframating_events listing mating events. Usually this is the output of
        `amajusmating.mating.mating_over_chains()`. At a minimum it should list
        columns 'mother', 'father', 'offspring' (number of offspring in each
        family).
    genotypes: faps.genotypeArray object
        Genotype data for plants in the mating pool including mothers and pollen
        donors.
    mu: float
        An estimate of the per-locus genotyping error rate.
        
    Returns
    =======
    A dictionary of FAPS genotypeArray objects for each maternal family.
    """
    # Pair of lists giving giving parents for each individual offspring
    noffspring = mating_events['offspring'].astype(int)
    dam_list  = [ [d]*n for d,n in zip(mating_events['mother'], noffspring) ]
    sire_list = [ [s]*n for s,n in zip(mating_events['father'], noffspring) ]
    # Flatten
    dam_list  = [x for y in dam_list  for x in y]
    sire_list = [x for y in sire_list for x in y]

    # Generate mating events, and generate arrays of progeny
    sim_progeny = fp.make_offspring(
        genotypes,
        dam_list = dam_list,
        sire_list= sire_list,
        mu = mu)
    # Add mutations to parental and simulated data.
    genotypes = genotypes.mutations(mu)
    sim_progeny = sim_progeny.mutations(mu)
    # Split progeny genotypes into maternal families
    sim_progeny = sim_progeny.split(by=sim_progeny.mothers)

    return sim_progeny

def simulate_dataset(data, model, mating_events, genotypes, mu, nloci):
    """
    Simulate a set of mating events, generate offspring genotypes, and calculate
    paternity arrays.

    Parameters
    ==========
    data: `faps_data` class object
        Data about the population.
    model: dict
        Dictionary of starting model parameters. Keys should be a subset of
        ['missing', 'shape', 'scale', 'mixture', 'assortment'] and values 
        floats giving initial values for those parameters, within appropriate
        boundaries.
    mating_events: pandas.DataFrame
        Dataframe listing mating events. Usually this is the output of
        `amajusmating.mating.mating_over_chains()`. At a minimum it should list
        columns 'mother', 'father', 'offspring' (number of offspring in each
        family).
        genotypes: faps.genotypeArray object
        Genotype data for plants in the mating pool including mothers and pollen
        donors.
    genotypes: faps.genotypeArray object
        Genotype data for plants in the mating pool including mothers and pollen
        donors.
    mu: float
        An estimate of the per-locus genotyping error rate.
    
    Returns
    =======
    Dictionary giving 
        'mating_events' : Dataframe of mating events giving mother, father and
            number of offspring.
        'genotypes' : Dict of faps.genotypeArray object for simulated offspring
            for each mating event
        'paternity_array' : Dict of faps.paternityArray objects for each 
            maternal family.
    """
    sim_data = {}
    # Table of simulated mating events giving the mother, father and family size
    sim_data['mating_events'] = pd.DataFrame({
        'mother' : mating_events['mother'],
        'father' : draw_sires(data, model, mating_events),# Draw candidate pollen donors for each mother.
        'offspring': mating_events['offspring']
    })

    # Simulate genotypes
    # Subset the correct number of loci
    locus_ix = np.sort(np.random.randint(genotypes.nloci, size=nloci))
    genotypes = genotypes.subset(loci=locus_ix)
    # Genotypes for the offspring, with point mutations
    sim_data['genotypes'] = simulate_genotypes(
        sim_data['mating_events'], genotypes, mu
    )

    # Probabilities of paternity for each offspring
    sim_data['paternity_array'] = fp.paternity_array(
        offspring = sim_data['genotypes'],
        mothers = { k: genotypes.subset(individuals = k) for k in sim_data['genotypes'].keys() }, # genotypes of the mothers
        males = genotypes,
        mu = mu,
        missing_parents = model['missing']
    )
    # Add dispersal probabilities as covariates
    for (p,s) in zip(sim_data['paternity_array'].keys(), data.covariates['dispersal']):
        sim_data['paternity_array'][p].add_covariate(s)

    # Add a vector of candidate-father names, including an entry for missing fathers.
    sim_data['candidates'] = np.array(data.candidates)
    sim_data['candidates'] = np.append(sim_data['candidates'], "nan")

    return sim_data

def remove_fathers(sim_data, prop_to_purge):
    """
    Remove a random subset of real fathers from the analysis.

    This samples real fathers without replacement and sets their probabilities
    of paternity for all offspring across all families to zero.

    Parameters
    ==========
    sim_data: dict
        The output of simulate_dataset
    prop_to_purge: float between zero and one
        Proportion of true fathers to be purged from the dataset.
        Note that this is the proportion of candidates that actually contributed
        to the progeny pool, not the proportion of overall candidates.

    Returns
    =======
    Nothing directly; the attribute for missing fathers is changed in each 
    paternityArray.
    """
    # Create a list of randomly selected fathers to remove from the data.
    unique_fathers = pd.unique(sim_data['mating_events']['father'])
    n_to_purge = np.round(
        len(unique_fathers) * prop_to_purge
    )
    fathers_to_purge = np.random.choice(unique_fathers, size = int(n_to_purge), replace=False)

    for k in sim_data['paternity_array'].keys():
        sim_data['paternity_array'][k].purge = fathers_to_purge # Set log likelihoods for purged candidates to -Inf

def mating_events_paternity_only(sim_data):
    """
    Table of mating events inferred from paternity only.

    This takes the most likely candidate father of each individual, including
    missing fathers, and returns a data frame of unique candidates and the 
    number of offspring belonging to each.

    Parameters
    ==========
    sim_data: dict
        The output of simulate_dataset

    Returns
    =======
    A dataframe giving the mother, inferred father and inferred number of offspring
    """
    inferred_mating_events = []
    for k in sim_data['paternity_array'].keys():
        max_ix = np.argmax(sim_data['paternity_array'][k].prob_array(), 1)
        top_candidates = sim_data['candidates'][max_ix]

        mating_table = pd.DataFrame({
            'mother'    : k,
            'father'    : np.unique(top_candidates),
            'offspring' : np.unique(top_candidates, return_counts = True)[1]
        })
        inferred_mating_events.append(mating_table)

    inferred_mating_events = pd.concat(inferred_mating_events)

    return inferred_mating_events


def accuracy_of_paternity(real_mating_events, inferred_mating_events):
    """
    Compare mating-event reconstruction by three methods for a single table of
    inferred mating events.    

    Parameters
    ==========
    real_mating_events: pd.DataFrame
        Table of correct mating events from simulate_dataset.
    inferred_mating_events: pd.DataFrame
        Table of inferred mating events.

    Returns
    =======
    A list of three floats giving:
        1. the probability that a mating event with >=1 offspring is correct
        2. the probability that a mating event with < offspring is correct
        3. the proportion of offspring assigned to a missing father.
    """

    pd.set_option('mode.chained_assignment', None)

    merge_inferred_and_real_events = pd.merge(
        real_mating_events, inferred_mating_events, 
        how='right', on=['mother', 'father'],
        suffixes = ['_real', '_inf']
        )
    
    # Get the probabilities that mating events with more or less than one offspring are real
    # Filter mating events with an inferred father
    inferred_mating_events = merge_inferred_and_real_events.loc[merge_inferred_and_real_events['father'] != 'nan']
    # # Boolean vector indicating if the inferred event was correct.
    inferred_mating_events['is_real'] = inferred_mating_events['offspring_real'].notna()
    # Probabilities that families of more or less than one are correct
    true_if_ge1 = inferred_mating_events.loc[inferred_mating_events['offspring_inf'] >=1, 'is_real'].mean()
    true_if_lt1 = inferred_mating_events.loc[inferred_mating_events['offspring_inf'] < 1, 'is_real'].mean()

    # Get the proportion of offspring that are assigned to a missing father. 
    n_unassigned = merge_inferred_and_real_events.loc[
        merge_inferred_and_real_events['father'] == 'nan', "offspring_inf"
        ].sum()
    p_unassigned = n_unassigned / merge_inferred_and_real_events['offspring_inf'].sum()
    
    return [true_if_ge1, true_if_lt1, p_unassigned]



def test_paternity_power(sim_data, prop_to_purge):
    """
    Assess the accuracy of paternity for an entire simulated dataset.

    This first purges a set of fathers, then infers mating events based on:
    1. paternity alone (taking the most-probable candidate of each progeny)
    2. paternity plus sibships; and paternity,
    3. sibships and phenotype information together.

    This then runs accuracy_of_paternity on each set of inferred mating events.

    Parameters
    ==========
    sim_data: dict
        The output of simulate_dataset
    prop_to_purge: float between zero and one
        Proportion of true fathers to be purged from the dataset.
        Note that this is the proportion of candidates that actually contributed
        to the progeny pool, not the proportion of overall candidates.
    
    Returns
    =======
    Dataframe with a row for the three inferred datasets:
        1. The probability that an inferred mating event is real given that there 
        were one or more offspring
        2. The probability that an inferred mating event is real given that there 
        were less than one offspring
        3. The proportion of offspring assigned to a missing father.
    """

    # Set the likelihood for a sample of fathers to zero
    remove_fathers(sim_data, prop_to_purge)

    # Cluster into sibships with and without covariate information.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore") # Suppress warnings about singleton sibships
        sibship_clusters_genetics_only   = fp.sibship_clustering(sim_data['paternity_array'], use_covariates=False)
        sibship_clusters_with_covariares = fp.sibship_clustering(sim_data['paternity_array'], use_covariates=True)

    # Summarise mating events
    # Return table giving the mother, inferred father, and number of offspring for 
    # 1.paternity only
    # 2. sibships + paternity
    # 3. paternity, sibships and covariate information
    inferred_mating_tables = {
        'paternity' : mating_events_paternity_only(sim_data),
        'sibships'  : fp.summarise_sires(sibship_clusters_genetics_only),
        'full'      : fp.summarise_sires(sibship_clusters_with_covariares)
    }

    # Calculate the probabilities that an inferred mating event is real, and the
    # proportion of offspring assigned to a missing father.
    # Returns three lists of length three, which as then flattened.
    sim_results = []
    for k in inferred_mating_tables.keys():
        sim_results.append(
            accuracy_of_paternity(sim_data['mating_events'], inferred_mating_tables[k])
        )
    sim_results = pd.DataFrame(
        np.array(sim_results),
        columns = ['true_ge1', 'true_lt1', 'missing']
        )
    sim_results.insert(0, 'data_type', inferred_mating_tables.keys())

    return sim_results