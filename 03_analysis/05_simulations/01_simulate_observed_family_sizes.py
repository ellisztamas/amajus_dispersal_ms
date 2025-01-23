"""
Simulate and analyse datasets using the observed parental genotypes using 
family sizes for each mother estimated for the data.

Returns a data frame giving:
1. The probability that an inferred mating event is correct, given that family
    size is >= 1.
2. The probability that an inferred mating event is correct, given that family
    size is < 1.
3. The proportion of offspring assigned to missing fathers.
"""

from amajusmating import simulations as sim
from tqdm import tqdm

np.random.seed(27)

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
outfile="03_analysis/05_simulations/simulate_observed_family_sizes.csv"
pd.DataFrame(
    {}, columns=['rep', 'scale', 'nloci', 'prop_purged', 'data_type', 'true_ge1', 'true_lt1', 'missing']
    ).\
        to_csv(outfile, index=False)

# Simulation parameters
nreps = 30
scale_values = [3, 30, 300]
nloci_values = [40, 60, 80]
q_values = [0.1, 0.3, 0.5]
total_simulations = nreps * len(scale_values) * len(nloci_values) * len(q_values)

sim_results = []
with tqdm(total=total_simulations) as pbar:
    for rep in range(nreps):
        for scale in scale_values:
            model={
                "missing" : 0.5,
                "mixture" : 0.9,
                "scale"   : scale,
                "shape"   : 1
            }
            for nloci in nloci_values:
                sim_data = sim.simulate_dataset(am_data, model, mating_events, adults, mu, nloci)
                for q in q_values:
                    pbar.update(1)
                    this_result = sim.test_paternity_power(sim_data, q)
                    this_result.insert(0, 'rep', rep)
                    this_result.insert(1, 'scale', scale)
                    this_result.insert(2, 'nloci', nloci)
                    this_result.insert(3, 'prop_purged', q)
                    this_result.to_csv(outfile, mode='a', index=False, header=False, float_format='%.3f')

