import numpy as np
import pandas as pd
import faps as fp

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

    return [x for y in sires for x in y]

def simulate_progeny(data, model, mating_events, genotypes, mu):
    """
    Simulate a genotypeArray of progeny.
    
    For each mother, this draws a sample of pollen donors in proportion to their
    distance from the mother, mates them, and simulates progeny genotypes.

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
    """
    # Vector of integers giving the number of offspring in each maternal family
    offspring_numbers = np.array(mating_events['offspring']).astype(int)
    
    # Draw candidate pollen donors for each mother.
    sire_list = draw_sires(data, model, mating_events)
    # Replicate each element to match the number of sires for each mother.
    sire_list = [ [s] * n for s,n in zip(sire_list,offspring_numbers) ]
    dam_list  = [ [m] * n for m,n in zip(mating_events['mother'], offspring_numbers) ]
    # Flatten the lists of maternal and paternal IDs
    sire_list = [x for y in sire_list for x in y]
    dam_list  = [x for y in dam_list for x in y]

    # simulate some progeny
    sim_progeny = fp.make_offspring(genotypes, dam_list=dam_list, sire_list=sire_list, mu = mu)
    
    return sim_progeny

def simulation_paternity(data, model, mating_events, genotypes, mu):
    """
    Get probabilities of paternity on simulated data.
    
    This simulates mating events between mothers and real candidate father
    genotypes, and generates families of full siblings of the same sizes as in 
    the real data sat. This then runs FAPS functions on those data and returns
    probabilities of paternity for true sires based on paternity alone, after
    clustering into sibships, and after clustering including information about
    dispersal.
    
    Paramating_eventsters
    ==========
    data: `faps_data` class object
        Data about the population.
    model: dict
        Dictionary of starting model paramating_eventsters. Keys should be a subset of
        ['missing', 'shape', 'scale', 'mixture', 'assortmating_eventsnt'] and values 
        floats giving initial values for those paramating_eventsters, within appropriate
        boundaries.
    mating_events: pandas.DataFramating_events
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
    A dataframe giving log probabilities of paternity for each offspring based on:
    1. Paternity only
    2. Sibship clustering
    3. Sibships plus covariates
    4. Probability that the father was missing when he was present
    5. Probability that the father was missing when he was missing
    """
    # Generate mating events, and generate arrays of progeny
    sim_progeny = simulate_progeny(
        data = data,
        model = model,
        mating_events= mating_events,
        genotypes = genotypes,
        mu = mu)
    
    # Add mutations to parental and simulated data.
    genotypes = genotypes.mutations(mu)
    sim_progeny = sim_progeny.mutations(mu)
    
    # Split progeny genotypes into maternal families
    sim_progeny = sim_progeny.split(by=sim_progeny.mothers)
    
    # Probabilities of paternity for each offspring
    patlik = fp.paternity_array(
        offspring = sim_progeny,
        mothers = {k: genotypes.subset(individuals=k) for k in sim_progeny.keys()}, 
        males = genotypes,
        mu = mu,
        missing_parents = model['missing']
    )
    # Add dispersal probabilities as covariates
    for (p,s) in zip(patlik.keys(), data.covariates['dispersal']):
        patlik[p].add_covariate(s)

    # Sibship clustering without covariates
    sc_raw = fp.sibship_clustering(patlik, use_covariates=False)
    # Sibship clustering with covariates
    sc_cov = fp.sibship_clustering(patlik, use_covariates=True)

    # Data.framating_events giving log probabilities of paternity for each offspring based on:
    # 1. Paternity only
    # 2. Sibship clustering
    # 3. Sibships plus covariates
    # 4. Probability that the father was missing when he was present
    output = {}
    for k in sim_progeny.keys():
        # Index positions of the true father
        ix = [np.arange(sim_progeny[k].size), sim_progeny[k].parent_index('father', genotypes.names)]
        
        output[k] = pd.DataFrame({
            'mother' : sim_progeny[k].mothers,
            'father' : sim_progeny[k].fathers,
            'paternity' : patlik[k].prob_array()[ ix[0], ix[1] ],
            'sibships'  : sc_raw[k].posterior_paternity_matrix()[ ix[0], ix[1] ],
            'sibs_covs' : sc_cov[k].posterior_paternity_matrix()[ ix[0], ix[1] ],
            'false_negative' :sc_cov[k].posterior_paternity_matrix()[:, -1]
        })

    output = pd.concat(output.values())
    
    return output
