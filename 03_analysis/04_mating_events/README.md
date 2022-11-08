Python scripts to infer mating events between mothers and individual pollen
donors. This is done for the output of the MCMC analysis in
`03_analysis/03_mcmc/01_mcmc_restrict_kurtosis/`. This takes each of the four
chains, discards the first 1500 iterations as burn-in, and infers mating events
based on dispersal parameters for each of the 1000 remaining iterations.

`get_mating_events()` calls the module
`02_library/python/amajusmating/mating.py`, which in turn is calling
`summarise_sires()` from FAPS using observed genotype and dispersal
data. For each maternal family, this pulls out a list of the plausible
fathers associated with each possible partition structure, and the posterior
probability of that structure. The probability of a mating event between the
mother and an individual candidate father is then the sum of those 
probabilities for each partition structure in which he appears.

`get_mating_events()` generates a file giving mating events between each mother
and plausible fathers for each iteration separately. This shows:

* ID of the mother
* ID of the father, 
* posterior support for that mating event,
* posterior mean number of offspring generated
* GPS and flower-colour information

`get_mating_events()` also generates a table summarising mating events for each
iteration.

`job_submission.sh` is a bash script that runs the mating script on each MCMC
output, and submits this as a job array to SLURM.