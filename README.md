# amajus_dispersal_ms
Code to support the manuscript "Joint estimation of paternity, sibships and pollen dispersal in a snapdragon hybrid zone"

##

## Folder overview

## Package dependencies

### Python

A conda environment file `environment.yml` is provided to make sure other users can use the same packages as were used to run the original analyses.
Assuming `conda` is installed on your machine (instructions [here](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)), install the environment with

```
conda env create -f environment.yml
```

Activate it before running analyses with
```
conda activate faps

```

### R

The cumulative pollen dispersal kernel is plotted with a custom function by 
[Nicolas Woloscko](https://github.com/NicolasWoloszko/stat_ecdf_weighted/blob/master/stat_ecdf_weighted.R)
to allow for weighted samples. This requires the CRAN package `spatstat`.