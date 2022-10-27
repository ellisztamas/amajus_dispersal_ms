import numpy as np
import pandas as pd
from tqdm import tqdm
import os

from amajusmating import simulations as sim
from amajusmating.mating import import_mcmc
# FAPS objects and distance matrices are generated in a separate script.
exec(open('03_analysis/01_data_formatting/setup_FAPS_GPS.py').read())

# Import the mating events we generated in a previous analysis.
mating_events = pd.read_csv("03_analysis/04_mating_events/output/mating_events_over_chains.csv")
# Trim mating events with probability < 0.9, and at least two offspring so we
# can test the effect of sibships clustering
mating_events = mating_events.\
loc[(mating_events['prob'] > 0.9) & 
    (mating_events['offspring'] > 1)]


# Parameters for posterior simulations of mating
np.random.seed(87)
burnin = 1500

# Import MCMC results with mating parameters at each iteration
input_dir = "03_analysis/03_mcmc/01_mcmc_restrict_kurtosis/output/"
mcmc = import_mcmc(input_dir, burnin=burnin)

# File to output the results.
output_file = "03_analysis/05_simulations/simulation_output.csv"
# We will append the file at each iteration, so first check whether there's a
# file there already and remove it.
if(os.path.isfile(output_file)):
    os.remove(output_file)
    print("Removing a previous output file.")

# Loop over steps in the MCMC chain and simulate mating events for each.
for i in tqdm(mcmc.index):
    model = mcmc.loc[i] # Get current dispersal model settings
    me = mating_events.loc[mating_events['iter'] == i] # Mating events for this iteration

    # Run the simulation
    this_sim = sim.simulation_paternity(am_data, model, me, adults, mu = 0.0013)
    this_sim.insert(loc=0, column='iter', value=i) # add a column giving the iteration label.
    
    # Write to disk
    with open(output_file, 'a') as f:
        this_sim.to_csv(f, mode='a', header=f.tell()==0, float_format='%.4f', index=False)