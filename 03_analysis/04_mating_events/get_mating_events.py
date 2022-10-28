#!/usr/bin/env python3

"""
Tom Ellis, 1st October 2021

Script to infer mating events for each iteration of an MCMC output. 
"""
import numpy as np
import os
import argparse

from amajusmating import mating


# Parameters
parser = argparse.ArgumentParser(description = 'Parse parameters to run the script')
parser.add_argument('-i', '--input', help = 'Path to input MCMC file. ', required = True)
args = parser.parse_args()

# FAPS objects and distance matrices are generated in a separate script.
exec(open('03_analysis/01_data_formatting/setup_FAPS_GPS.py').read())

# Parameters for posterior simulations of mating
np.random.seed(87)
burnin = 1500

# Input and output data
input_dir = args.input
output_dir = os.path.dirname(os.path.abspath(__file__)) + '/output/' + os.path.basename(args.input) + '/'
os.makedirs(output_dir, exist_ok=True)

print("\nInferring mating events from each iteration of the MCMC chain, discarding the first {} iterations as burn-in.\n\n".format(burnin))
mating.mating_over_chains(
    data = am_data,
    input_dir = input_dir + '/output/',
    output_dir = output_dir,
    burnin = burnin
    )
