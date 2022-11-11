This folder contains the following data files:

* `GPS_positions.csv`: Raw GPS data for all flowering plants sampled in 2012
* `offspring_2012_genotypes.csv`: SNP data for the offspring
* `parents_2012_genotypes.csv`: SNP data for flowering plants sampled in 2012, 
    including the maternal plants
* `paternity_array.csv`: A matrix of probabilities of paternity for each
    candidate father on each offspring. This is created when you run the script
    `03_analysis/01_data_formatting/setup_FAPS_GPS.py` and is not included in 
    the GitHub repo. This is just to save time when calling the script later
    for other analyses.
* `processed_GPS_positions.csv`: The same data as `GPS_positions.csv` but
    formatted by `03_analysis/01_data_formatting/setup_FAPS_GPS.py`. It is not
    included in the GitHub repo, but is created by that script. This is
    here so we can use GPS data for plotting in R. 
* `rosea_sulfurea.csv`: Flower-colour phenotype data for each plant sampled in 
    2012.