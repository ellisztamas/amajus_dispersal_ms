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

def simulate_dataset(data, model, mating_events, genotypes, mu):
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

    return sim_data

def simulate_mating_events(sim_data, q_real=0, q_param=None):
    """
    Infer mating events from simulated data.

    This takes a dictionary of correct mating events and the corresponding
    paternity array, removes a sample of fathers at random, and infers mating
    events from the data.   

    Parameters
    ==========
    sim_data: dict
        Dictionary giving a set of correct mating events and the corresponding
        paternity array. This is the output of `simulate_data()`.
    q_real: float
        True proportion of missing fathers.
    q_param: float
        Input value for the proportion of missing fathers used as a prior for 
        paternity inference. Defaults to q_real.
    
    Returns
    =======
    Dataframe listing all correct and inferred mating events, with columns:
        'mother' : Name of the mother
        'father' : Name of the father.
        'offspring_real' : Correct number of offspring. For incorrectly inferred
            mating events this is NaN.
        'prob' : Posterior probability of the inferred mating event.
        'offspring_inf' : Number of inferred offspring. NaN for correct mating 
            events that were not detected.
        'q_real' : True proportion of missing fathers.
        'q_param': Input value for the proportion of missing fathers used as a prior for 
            paternity inference.
    """
    if q_param is None:
        q_param = q_real
        
    # Create a list of randomly selected fathers to remove from the data.
    n_to_purge = int(np.round(len(sim_data['mating_events']['father']) * q_real))
    fathers_to_purge = np.random.choice(sim_data['mating_events']['father'], size =  n_to_purge, replace = False)

    for k in sim_data['paternity_array'].keys():
        sim_data['paternity_array'][k].missing_parents = q_param # Set the prior proportion of missing fathers
        sim_data['paternity_array'][k].purge = fathers_to_purge # Set log likelihoods for purged candidates to -Inf

    # Cluster into sibships, and get a list of inferred mating events
    with warnings.catch_warnings():
        warnings.simplefilter("ignore") # Suppress warnings about singleton sibships
        sibship_clusters = fp.sibship_clustering(sim_data['paternity_array'], use_covariates=False)
        inferred_mating_events = fp.summarise_sires(sibship_clusters)
    # Merge the real and inferred mating events
    # Returns a dataframe for all mating events including real and incorrect false ones
    output= sim_data['mating_events'].merge(
        inferred_mating_events, how="outer", on = ('mother', 'father'), suffixes = ('_real', '_inf')
        )
    output = output.drop('log_prob', axis=1)
    output['q_real'] = q_real
    output['q_param'] = q_param


    return output
    


# def simulate_missing_fathers(data, progeny, paternity_array, q):
    # """
    # Remove a subset of real fathers from a paternityArray object and see what 
    # happens to the paternity of the offspring.

    # This relies on faps.summarise_paternity() to summarise the top candidate for
    # each progeny and counts how often fathers of each offspring are 
    # correctly/incorrectly assigned to be missing.

    # Parameters
    # ==========
    # data: `faps_data` class object
    #     Data about the population.
    # progeny: faps.genotypeArray object
    #     Genotype data on a set of progeny
    # paternity_array: faps.paternityArray object
    #     Inferred paternity on the full data set
    # q: float between zero and one
    #     Proportion of real fathers to remove from the data by adding to the 
    #     `purge` attribute of the paternityArray.

    # Returns
    # =======
    # Dictionary giving:
    #     'q'                     : Proportion of fathers removed
    #     'total_present':        : Number of offspring whose father was present
    #     'true_positive'         : Probability that a true father is the most likely candidate, assuming he is sampled
    #     'false_negative'        : Probability that a father is most likely missing, even if he is sampled
    #     'total_orphans'         : Number of progeny whose father was removed
    #     'identified_as_missing' : Number of orphaned progeny whose father was identified as missing
    #     'false_singleton'       : Number of orphaned progeny assigned to a singleton family
    #     'prob_false_paternity'  : Mean posterior probability of paternity for falsely assigned fathers
    #     'median_distance_real'  : Weighted-median dispersal distance in the complete dataset
    #     'median_distance_inferred' : Weighted-median dispersal distance after subsampling fathers

    # """
    # # Purge some of the real fathers and cluster into sibships again
    # to_purge = [ k.fathers for k in progeny.values() ]
    # to_purge = [ x for y in to_purge for x in y ]
    # to_purge = np.unique(to_purge)
    # to_purge = np.random.choice(to_purge, size = int(np.round(to_purge.shape[0] * q)) )

    # # Dataframe giving mother, father, probability of paternity and offspring number for the purged data
    # for k in paternity_array.keys():
    #     paternity_array[k].purge = to_purge
    # sc_purged = fp.sibship_clustering(paternity_array, use_covariates=True)
    # me_purged = fp.summarise_paternity(sc_purged)
    # # Data frame giving progeny and the real father
    # me_real = {}
    # for k in progeny.keys():
    #     me_real[k] = pd.DataFrame({
    #         'progeny' : progeny[k].names,
    #         'real_father' : progeny[k].fathers
    #         })
    # me_real = pd.concat(me_real)

    # # Merge real and inferred fathers
    # me_combined = pd.merge(me_real, me_purged, on="progeny")
    # # subset those offspring whose father was/was not missing
    # me_purged = me_combined[me_combined['real_father'].isin(to_purge)]
    # me_not_purged = me_combined[~me_combined['real_father'].isin(to_purge)]

    # # Number of progeny whose father was sampled or removed
    # total_present = me_not_purged.shape[0]
    # total_orphans = me_purged.shape[0]
    # # Frequency of missing fathers inferred to be missing
    # orphans_identified_as_missing = np.sum( me_purged['candidate_1'] == "missing"  )
    # # Probability that true fathers are identified, given that they are sampled
    # true_positive = np.mean( me_not_purged['real_father'] == me_not_purged['candidate_1'] )
    # # Probability that true fathers are inferred to be missing
    # false_negative = np.mean( me_not_purged['candidate_1'] == "missing"  )  
    # # Probability that when the father is missing, offspring are falsely assigned to a family of size=1.
    # orphans_in_singleton_family = np.sum(
    #     me_purged[me_purged['candidate_1'] != "missing"].groupby(['mother', 'candidate_1']).size() == 1
    # )
    # # Mean posterior probability for false families
    # prob_false_paternity = np.mean(np.exp(me_purged[me_purged['candidate_1'] != "missing"]['logprob_1']) )

    # # Calculate weighted median distances for full and purged data
    # # Create a data frame of distances we can join on.
    # distance_df = pd.DataFrame({
    #     'mother'   : np.repeat(list(data.mothers), data.n_candidates),
    #     'father'   : np.tile(data.candidates, len(data.mothers)),
    #     'distance' : data.distances.flatten()
    # })
    # # Merge siring events with distances and phenotypes
    # median_distance_real = me_combined.\
    #     merge(distance_df, how="left", left_on=['mother', 'real_father'], right_on = ['mother', 'father']).\
    #         drop_duplicates(subset=['mother', 'real_father'])
    
    # median_distance_purged = me_purged.\
    #     merge(distance_df, how="left", left_on=['mother', 'candidate_1'], right_on = ['mother', 'father']).\
    #     drop_duplicates(subset=['mother', 'real_father'])
    # median_distance_purged = median_distance_purged.loc[~median_distance_purged['distance'].isna()]

    # median_distance_inferred = me_combined.\
    #     merge(distance_df, how="left", left_on=['mother', 'candidate_1'], right_on = ['mother', 'father']).\
    #     drop_duplicates(subset=['mother', 'real_father'])
    
    # # Calculate median dispersal distance weighted by posterior support for each mating event.
    # median_distance_real     = wquantiles.median(median_distance_real['distance'],   np.exp(median_distance_real['logprob_1']))
    # median_distance_inferred = wquantiles.median(median_distance_inferred['distance'], np.exp(median_distance_inferred['logprob_1']))
    # return {
    #     'q'                     : q,
    #     'total_present'         : total_present,
    #     'true_positive'         : true_positive,
    #     'false_negative'        : false_negative,
    #     'total_orphans'         : total_orphans,
    #     'identified_as_missing' : orphans_identified_as_missing,
    #     'false_singleton'       : orphans_in_singleton_family,
    #     'prob_false_paternity'  : prob_false_paternity,
    #     'median_distance_real'  : median_distance_real,
    #     'median_distance_inferred' : median_distance_inferred
    # }