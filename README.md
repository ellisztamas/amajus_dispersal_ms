# amajus_dispersal_ms
Code to support the manuscript "Joint estimation of paternity, sibships and pollen dispersal in a snapdragon hybrid zone"

## Folder overview

The folder is separated into folder with names that should be self-explanatory:

1. Data files, including a README explaining what they are.
2. Library code, including a python package to do the heavy lifting of running
    the MCMC, inferring mating events and running simulations. There are also 
    some R functions to do odd jobs
3. Analysis folders. These have their own README files, except for the one for 
    prior simulations, which has a RMarkdown file to explain what it does.
4. Code to create the figures. Actual figure files save to the manuscript folder.
5. Files to create the manuscript.

## To recreate the analyses

- All scripts are run fom the folder root.
- Assuming you have Python installed, install the conda environment (see below)
- Analyses are run in the order they appear in the folder.
- For most analyses I committed the output files to this repo. The exception are
the output files for the inference of mating events, because these were massive.
You will have to run these yourself to run any downstream analyses (e.g. 
simulations) or plots - there's a job submission script to do this.

## Package dependencies

### Python

A conda environment file `environment.yml` is provided to make sure other users 
can use the same packages as were used to run the original analyses.
Assuming `conda` is installed on your machine 
(instructions [here](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)),
install the environment with

```
conda env create -f environment.yml
```

Activate it before running analyses with
```
conda activate faps

```

### R

The file `session_info` gives the full list of package versions used, but in brief:

- The `tidverse` collection
- `cowplot` and the related `ggarrange`
- The cumulative pollen dispersal kernel is plotted with a custom function by 
[Nicolas Woloscko](https://github.com/NicolasWoloszko/stat_ecdf_weighted/blob/master/stat_ecdf_weighted.R)
to allow for weighted samples. This requires the CRAN package `spatstat`.
- `rmarkdown` and `knitr` to create the supplementary material.