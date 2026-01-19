#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 11:21:01 2024

@author: nrakotondrainibe
"""
import matplotlib.pyplot as plt
import numpy as np
import os
from astropy.io import ascii
from astropy.table import Table, vstack
import pandas as pd
import importlib.util
spec = importlib.util.find_spec("pyGRBz")
from pyGRBz.estimation import stats
if spec and spec.origin:
    pathgrbz = os.path.dirname(spec.origin)
    path = os.getcwd()
else:
    print("Warning: Could not find module 'pyGRBz'. Using current working directory.")
from pyGRBz.pyGRBz import GRB_photoZ    

grb_name = input("GRB name: ").replace(" ", "").upper() #format GRB050904
grb_name = [grb_name]

#grb_name = ["GRB050904"]
mode = input("SED or LC? ").lower() # Can estimate photo_z for: "MutlipleTargets, SED, LC"

# input data and output
if mode=="sed":
    input_dir='/zmodule/data/sed/'
elif mode=="lc":
    input_dir='/zmodule/data/lc/'
else:
    sys.exit("Wrong input")

#output_dir=path+'/results/%s/'%mode
output_dir='/zmodule/results/%s/'%mode


# Load module
photoz = GRB_photoZ(
    output_dir=output_dir,
    thres_err=0.02, # if flux_err/flux < thres_err then set flux_err = thres_err*flux
    wvl_step=50, # angstroms
    wvl_step_X=10 # angstroms
)

###############################################################################################################
# Load as many targets as you want. It can be a mix of SEDs and light curves
photoz.load_data(
    data_dir=input_dir,           
    data_name=grb_name
)

###############################################################################################################
#Format data in order to apply galactic extinction if MW not corrected and calculates the flux in Jansky to each observation, 
#E(B-V) values from Schlegel, Finkbeiner & Davis (1998) directly from https://irsa.ipac.caltech.edu/applications/DUST/, can use the calibration of SFD11.
#Deredenning using the Cardelli 1989 parametrization
photoz.formatting()


###############################################################################################################
# Extract the SED at a given time.
# First the data are fitted either with a single power law (SPL) or a broken power law (BPL): model = "SPL","BPL"
# Secondly the time at which to extract the SED can be either 'fixed' (needs to give through time_SED in seconds) or 
# computed to be the time at which the flux is maximum in the reddest band ('ReddestBand'):
# method='ReddestBand' or method='fixed',time_SED = ...

# In case the input data is already a SED. This function has to run in order to have the right formatting for the follwing computations
photoz.extract_sed(model='SPL',method='ReddestBand')


###############################################################################################################
# Create flat priors for the photo_z estimation
priors=dict(z=[0,11],Av=[0,10],beta=[0,3],norm=[0,10])#NHX if with X-data

# Run the MCMC algorithm.
# Select the extinction law to used: 'smc', 'lmc', 'mw', 'nodust'
# Nthreads: number of threads to use in case of parallelization
# nwalkers: number of walkers
# Nsteps1: number of steps for the first burn-in phase
# Nsteps2: number of steps for the second burn-in phase
# Nsteps3: number of steps for the production run
# Select to add dust, gas in host and our Galaxy
# Select IGM transmission method: 'Madau99' or 'Meiksin06'
# adapt_z: If adapt_z is True try to reduce parameter space for the redshift based on non detection in blue bands
laws = ['smc','nodust']
#laws = ['smc', 'lmc', 'mw', 'nodust','sne']
for law in laws:
    photoz.fit(
        ext_law=law,
        Nthreads=6,
        nwalkers=50,
#        Nsteps1=500,
#        Nsteps2=1000,
        Nsteps1=0,
        Nsteps2=1000,
        nburn=300,
#        nburn=300,
        Host_dust=True,
        Host_gas=False,
        igm_att='Meiksin',
        clean_data=False,
        priors=priors,
        adapt_z=True)

#output_dir=path+'/results/%s/'%mode
output_dir='/zmodule/results/%s/'%mode
###############################################################################################################    
#Compare all fits statistically, in term of chi-square and BIC
user_input = input("Do you want to statistically compare the results of your different fits? (yes/no)").replace(" ", "").lower()
if user_input.lower() in ["yes", "y"]:
	stats(path+output_dir+grb_name[0]+"/",ext_laws=laws,lim_bic=2)	
elif (user_input.lower() in ["no", "n"]) or len(laws)<2:
    print("You have only one model anyway...")
else:
    print("Invalid input but you have only one model anyway...")
