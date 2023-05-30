"""
Tom Ellis, October 2019.

This script imports and prepares data on
1. Genotype data of parents and offspring
2. GPS data
3. Flower colour phenotypes

It then creates and saves a paternity array object for genotype data.
If this file exists when this script is called again, the saved version is imported
Finally genotype and paternity arrays are split by maternal family.

It returns a faps_data object with all the relevant information.
"""

import numpy as np
import faps as fp
import pandas as pd
from os.path import isfile
# Import local modules
from amajusmating import faps_data

print(
    "\nRunning the script to format genotype, flower-colour and GPS information using FAPS version {}.\n".format(fp.__version__)
    )

# GENOTYPE DATA
# Genotyping error rate.
mu=0.0001
# Import genpotype data
print("Importing genotype data.")
progeny = fp.read_genotypes('01_data/offspring_2012_genotypes.csv', mothers_col=1, genotype_col=2)
adults  = fp.read_genotypes('01_data/parents_2012_genotypes.csv')
# SNP Data cleaning.
print("Filtering loci and individuals with poor marker data.")
# remove individuals with >7.5% missing data
md = 0.075
progeny = progeny.subset(individuals= progeny.missing_data('individual') < md)
adults  = adults.subset(individuals= adults.missing_data('individual') < md)
# remove loci with more than 10% missing data
r = 0.1
lx = (progeny.missing_data('marker') < r) * (adults.missing_data('marker') < r)
progeny = progeny.subset(loci=lx)
adults  = adults.subset(loci=lx)
# remove loci with heterozygosity less than 15%
h = 0.15
progeny = progeny.subset(loci= (adults.heterozygosity('marker') > h))
adults  = adults.subset( loci= (adults.heterozygosity('marker') > h))

# PATERNITY ARRAYS
# If this script is run for the first time, the paternityArray is saved for later.
# If the file exists already, this is imported, which is much faster.
# Otherwise a new paternityArray object is created.
if isfile("01_data/paternity_array.csv"):
    print(
        """Importing the matrix of probabilities of paternity. 
To create this file from scratch, manually delete `01_data/paternity_array.csv` and run this script again.""")
    # genotypeArray giving the mother of each maternal family.
    mothers = adults.subset(individuals=np.unique(progeny.mothers))
     # Split into maternal families.
    # Import saved paternityArray.
    patlik = fp.read_paternity_array("01_data/paternity_array.csv")
else:
    print("Creating and saving a matrix of probabilities of paternity...")
    # A single genotypeArray giving the mother of each of 984 offspring individuals.
    mothers = adults.subset(individuals=progeny.mothers)
    # Create the paternity array and save for later.
    patlik = fp.paternity_array(
        progeny, mothers, adults, mu = mu
        )
    patlik.write("01_data/paternity_array.csv")

# SPLIT BY MATERNAL FAMILYs
# genotypeArray giving the mother of each maternal family.
print("Splitting the paternityArray object into maternal families.")
mothers = adults.subset(individuals=np.unique(progeny.mothers))
# Split objects up by maternal family.
mothers = mothers.split(np.unique(progeny.mothers))
patlik  = patlik.split(by=progeny.mothers)
progeny = progeny.split(by=progeny.mothers)

# FILTER OUT FATHERS WITH OPPOSING GENOTYPES
# For each maternal family, create a matrix giving the number of loci at which 
# a candidate and offspring have opposing homozygous genotypes.
for k in patlik.keys():
    patlik[k].clashes = fp.incompatibilities(adults, progeny[k])
    patlik[k].max_clashes=2 # Define a maximum number of homozygous incompatibilities are allowed for each trio

# GPS DATA
print("Incorportating GPS and flower-colour data.")
# Import GPS data
gps = pd.read_csv("01_data/GPS_positions.csv", index_col=0)
gps = gps.loc[adults.names] # reorder to match cleaned SNP data
gps = gps[['Easting','Northing']] # remove the column for altitude. We don't need that here.
# Save GPS positions so they can be used in R.
gps.\
    assign(
        is_mother = gps.index.isin(list(mothers.keys())) # Add a column stating whether the plant is also a mother
    ).\
        to_csv("01_data/processed_GPS_positions.csv", index = True) 

# FLOWER COLOUR
# Import flower colour data for the population
ros_sulf = pd.read_csv('01_data/rosea_sulfurea.csv', index_col='id')
ros_sulf = ros_sulf.loc[adults.names] # ensure data are in the same order as for the genotype data
# Simplify flower colours to yellow, full red or hybrid
ros_sulf['simple_colour'] = 'unkown'
ros_sulf.loc[ros_sulf['flower_colour'].isin(['FR',"Ye"]), 'simple_colour'] = ros_sulf['flower_colour']
ros_sulf.loc[ros_sulf['flower_colour'].isin(['FO',"WR", "WO", "Wh"]), 'simple_colour'] = 'hybrid'

# Create a class object to hold the data.
print("Creating a `faps_data` object to hold genotype and phenotype information.\n")
am_data = faps_data(
    paternity=patlik,
    gps = gps,
    flower_colours = ros_sulf[['rosea','sulfurea', 'simple_colour']],
    params = {}
    )

# A table of maternal family sizes for each mother.
# Save it now so we can easily import in R later
pd.DataFrame({
    'mother' : am_data.mothers,
    'n_offspring': [len(x.offspring) for x in am_data.paternity.values()]
}).to_csv("01_data/maternal_family_sizes.csv", index = False)


# Remove variables we don't need anymore.
del(md, r, h, ros_sulf, patlik, lx)

print("Data formatting script completed.\n\n\n")