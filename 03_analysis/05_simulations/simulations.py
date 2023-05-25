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
burnin = 1500

# Import MCMC results with mating parameters at each iteration
input_dir = "03_analysis/03_mcmc/01_mcmc_main/output/"
mcmc = import_mcmc(input_dir, burnin=burnin)
# Just use 200 random iterations from the MCMC
ix = np.sort( np.random.random_integers(0,mcmc.shape[0], 200) )
mcmc = mcmc.iloc[ix, :]

# File to output the results.
output_dir = "03_analysis/05_simulations"
# We will append the file at each iteration, so first check whether there's a
# file there already and remove it.
if(os.path.isfile(output_dir + "/compare_paternity_accuracy.csv")):
    os.remove(output_dir + "/compare_paternity_accuracy.csv")
    print("Removing output file: " + output_dir + "/compare_paternity_accuracy.csv")
if(os.path.isfile(output_dir + "/missing_fathers.csv")):
    os.remove(output_dir + "/missing_fathers.csv")
    print("Removing output file: " + output_dir + "/missing_fathers.csv")


# Loop over steps in the MCMC chain and simulate mating events for each.
for i in tqdm(mcmc.index):
    model = mcmc.loc[i] # Get current dispersal model settings
    me = mating_events.loc[mating_events['iter'] == i] # Mating events for this iteration

    # Simulate progeny and calculate paternity likelihoods
    sim_progeny, patlik = sim.simulation_paternity(am_data, model, me, adults, mu = 0.0001)
        
    # Get the accuracy of paternity based on paternity only, paternity plus sibships, and the joint model
    accuracy = sim.compare_paternity_accuracy(sim_progeny, patlik, adults)
    accuracy.insert(loc=0, column='iter', value=i) # add a column giving the iteration label.
    # Write to disk
    with open(output_dir + "/compare_paternity_accuracy.csv", 'a') as f:
        accuracy.to_csv(f, mode='a', header=f.tell()==0, float_format='%.4f', index=False)
    
    # Remove a proportion of the true fathers and see how this affects paternity assignments.
    # See amajus.sims.missing_fathers for details of the output.
    missing_dads = [sim.simulate_missing_fathers(am_data, sim_progeny, patlik, q) for q in [0.1, 0.2, 0.3, 0.4, 0.5] ]
    missing_dads = pd.DataFrame(missing_dads)
    missing_dads.insert(loc=0, column='iter', value=i) # add a column giving the iteration label.
    # Write to disk
    with open(output_dir + "/missing_fathers.csv", 'a') as f:
        accuracy.to_csv(f, mode='a', header=f.tell()==0, float_format='%.4f', index=False)
    
