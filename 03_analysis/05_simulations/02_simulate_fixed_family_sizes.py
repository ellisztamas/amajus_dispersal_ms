"""
Simulate and analyse datasets using the observed parental genotypes using 
a fixed number of offspring per full-sib family.

Returns a data frame giving:
1. The probability that an inferred mating event is correct, given that family
    size is >= 1.
2. The probability that an inferred mating event is correct, given that family
    size is < 1.
3. The proportion of offspring assigned to missing fathers.
"""

from amajusmating import simulations as sim
from time import time
import numpy as np
import multiprocessing
import os


np.random.seed(23)

# Import and format raw data
exec(open('03_analysis/01_data_formatting/setup_FAPS_GPS.py').read())

# Import the mating events we generated in a previous analysis.
mating_events = pd.read_csv("03_analysis/04_mating_events/output/01_mcmc_main/mating_events_over_chains.csv")
# Get posterior mean posterior probabilities and sizes for each sibship
mating_events = mating_events.\
    groupby(['mother','father']).mean().\
    reset_index()
# Trim mating events with low support
mating_events = mating_events.loc[(mating_events['prob'] >= 1)]

# Empty CSV file to store the output
mating_file="03_analysis/05_simulations/simulate_fixed_family_sizes_mating.csv"
pd.DataFrame(
    {}, columns=[
        'rep', 'scale', 'nloci', 'prop_purged', 'family_size', 'data_type', 
        'arith_prob_correct', 'weighted_prob_correct', 
        'true_if_pp_is1', 'true_if_pp_lt1', 
        'true_if_pp_is99', 'true_if_pp_lt99',
        'true_if_offs_ge1', 'true_if_offs_lt1',
        'missing']
    ).\
        to_csv(mating_file, index=False)

paternity_file = "03_analysis/05_simulations/simulate_fixed_family_sizes_paternity.csv"
pd.DataFrame({}, columns = [
    'rep', 'scale', 'nloci', 'prop_purged', 'family_size', 'data_type', 
    'true_pos', "false_pos", "true_neg", "false_neg"
]). to_csv(paternity_file, index=False)

# Simulation parameters
nreps = 30
offs_values = [1,2,3]
scale_values = [3,30,300]
nloci_values = [40, 50, 60]
q_values = [0.1, 0.3, 0.5]
total_simulations = nreps * len(offs_values) * len(scale_values) * len(nloci_values) * len(q_values)

t0 = time()

def run_iteration(rep):
    print("Starting replicate simulation {}.".format(rep))
    t1 = time()
    for noffs in offs_values:
        for scale in scale_values:
            model={
                "missing" : 0.5,
                "mixture" : 0.9,
                "scale"   : scale,
                "shape"   : 1
            }
            for nloci in nloci_values:
                # Ensure each parallel process starts with a different seed
                np.random.seed((os.getpid() * int(time())) % 123456789)

                sim_data = sim.simulate_dataset(
                    am_data, model, mating_events, adults, mu, nloci, noffspring=noffs
                    )
                
                for q in q_values:
                    # Calculate accuracies
                    sim_result = sim.test_paternity_power(sim_data, q)
                    for v in sim_result.values():            
                        v.insert(0, 'rep', rep)
                        v.insert(1, 'scale', scale)
                        v.insert(2, 'nloci', nloci)
                        v.insert(3, 'prop_purged', q)
                        v.insert(4, 'family_size', noffs)
                    # Save to disk
                    sim_result['mating'].to_csv(   mating_file,    mode='a', index=False, header=False, float_format='%.3f')
                    sim_result['paternity'].to_csv(paternity_file, mode='a', index=False, header=False, float_format='%.3f')
                    
    t2 = time()
    time_for_this_rep = np.round((t2-t1)/60, 2)
    time_so_far = np.round((t2-t0)/3600, 2)
    print("Replicate {} completed in {} seconds. {} hours have elpased so far.".format(rep, time_for_this_rep, time_so_far))

if __name__=='__main__':
    pool = multiprocessing.Pool()
    results = pool.map(run_iteration, range(nreps))
    pool.close()
    pool.join()
