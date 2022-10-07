This folder contains scripts to run Metropolis-Hastings MCMC for dispersal, 
sibships, and paternity.

There are five folders corresponding to analyses with different priors. Unless
otherwise stated, analyses used a prior proportion of missing fathers of 0.32.
For more details on these prior settings, see
`03_analysis/03_prior_simulations/prior_simulations.Rmd`,but here is a summary:

1. Priors that are moderately skeptical about leptokurtic dispersal
2. Priors favouring leptokurtosis at short scale, but skeptical about long-range
    migration.
3. Priors that strongly penalise leptokurtosic dispersal.
4. The same as 1, but with the proportion of missing fathers set to 0.22
5. The same as 1, but with the proportion of missing fathers set to 0.42

Each folder contains a file called `run_MCMC.py` to run the analysis and save
the output. This is mostly calling library code from
 `02_library/python/amajusmating/mcmc.py` - look there for the nuts and bolts.
 The `run_MCMC.py` scripts contain (hard-coded) dictionaries giving:

1. Prior settings
2. Sigma values for the proposal distribition. These are standard deviations of
    a normal distribution used to peturb each variable at each iteration.
3. Initial parameter values to initialise the chains

Each script runs four independent chains with 4000 iterations. There are no
burn-in iterations here - these are pruned later.

If you have access to a computing cluster which uses SLURM, there is also a
script `job_submission.sh` which can be used to submit all the `run_MCMC.py`
files as a batch job. I have tested this on the GMI cluster, and it worked, but 
it may need tweaking to get it working on other systems.

Note that this MCMC is not very fast, because you have to perform the 
Monte-Carlo simulations for sibship clustering at each iteration. On both my
laptop and the GMI computing cluster this took 3-4 seconds per iteration. I am
sure there is scope to speed this up should someone with a stronger computer
science background than wish to do so.

Tom Ellis, 7th October 2022