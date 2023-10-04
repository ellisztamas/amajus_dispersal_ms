import numpy as np
import pandas as pd
from tqdm import tqdm
import os

from amajusmating import simulations as sim
from amajusmating.mating import import_mcmc
# FAPS objects and distance matrices are generated in a separate script.
exec(open('03_analysis/01_data_formatting/setup_FAPS_GPS.py').read())

# Import the mating events we generated in a previous analysis.
mating_events = pd.read_csv("03_analysis/04_mating_events/output/01_mcmc_main/mating_events_over_chains.csv")
# Trim mating events with probability < 0.9, and at least two offspring so we
# can test the effect of sibships clustering
mating_events = mating_events.\
loc[(mating_events['prob'] > 0.9) & 
    (mating_events['offspring'] > 1)]

# Parameters for posterior simulations of mating
np.random.seed(87)
burnin = 500

# Import MCMC results with mating parameters at each iteration
input_dir = "03_analysis/03_mcmc/01_mcmc_main/output/"
mcmc = import_mcmc(input_dir, burnin=burnin)
# Just use 200 random iterations from the MCMC
ix = np.sort( np.random.randint(0,mcmc.shape[0], 200) )
mcmc = mcmc.iloc[ix, :].reset_index()

# File to output the results.
output_dir = "03_analysis/05_simulations"
# The simulation appends the file at each iteration, so first check whether there's a
# file there already and remove it.
if(os.path.isfile(output_dir + "/simulated_mating_events.csv")):
    os.remove(output_dir + "/simulated_mating_events.csv")
    print("Removing output file: " + output_dir + "/simulated_mating_events.csv")

# Loop over steps in the MCMC chain and simulate mating events for each.
for i in tqdm(mcmc.index):
    model = mcmc.loc[i] # Get current dispersal model settings
    me = mating_events.loc[mating_events['iter'] == i] # Mating events for this iteration

    # Simulate a dictionary of mating events, offspring genotypesand paternity likelihoods
    sim_data = sim.simulate_dataset(am_data, model, me, adults, mu)

    # Infer mating events for different values of proportion of missing fathers (q)
    sim_mating_events = [ sim.simulate_mating_events(sim_data, q) for q in [0.1, 0.2, 0.3, 0.4, 0.5] ]
    sim_mating_events = pd.DataFrame(sim_mating_events)

    # Write to disk
    with open(output_dir + "/simulated_mating_events.csv", 'a') as f:
        sim_mating_events.to_csv(f, mode='a', header=f.tell()==0, float_format='%.4f', index=False)
    
