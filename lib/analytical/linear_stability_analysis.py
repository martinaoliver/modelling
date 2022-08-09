#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 12:33:13 2020

@author: mo2016
"""
# - This file calls the findsteadystates_functions.py and dispersionrelation_functions.py to carry
# linear stability analysis on a desired parameter set. The parameter set can be inputed either as a
# dictionary (single parameter set) or as a dataframe (multiple parameter sets).


from analytical.findsteadystates_functions import findsteadystates
from analytical.dispersionrelation_functions import dispersionrelation
import pandas as pd
import numpy as np

#Turing analysis carried out on a dataframe. The input is a df with every parameter set.
def big_turing_analysis_df(df,circuit_n,n_species,top_dispersion=5000,print_parID=False):
    len_df = len(df) #lenght of dataframe (number of parameter sets to analyse)
    output_df = pd.DataFrame(data=None, columns=df.columns)
    # par_dict['ss_n'],par_dict['ss_list'],par_dict['ss_class'],par_dict['system_class'],par_dict['maxeig'],par_dict['new_index'] =[np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]

    for parID in df.index:
        if print_parID == True:
            print(parID)
        
        par_dict = df.loc[parID].to_dict() #converts a dataframe row into a dictionary outputing a dictionary for a specific parameter set
        steadystatelist, number_steadystates = findsteadystates(par_dict,circuit_n,n_species, n_initial_conditions = 100) #input a dictionary with the parameters and returns (1) a list with the steady states and (2) the number of steady states.
        if number_steadystates > 0:
            for ss_n in range(number_steadystates): #perform linear stability analysis on all steady states found
                steadystate_values_ss_n = steadystatelist[ss_n]
                ss_class, system_class, eigenvalues, maxeig= dispersionrelation(par_dict,steadystate_values_ss_n, circuit_n,top_dispersion)
                # maxeig, pattern_class, eigenvalues ,oscillations, eigenvalsteadystate= dispersionrelation(par_dict,steadystate_values_ss_n, circuit_n)
                # par_dict['ss_n'],par_dict['ss_list'],par_dict['class'],par_dict['maxeig'],par_dict['oscillations'],par_dict['k0_eig'],par_dict['new_index'] = number_steadystates,steadystate_values_ss_n,pattern_class,maxeig,oscillations,eigenvalsteadystate,[parID,ss_n]
                par_dict['ss_n'],par_dict['ss_list'],par_dict['ss_class'],par_dict['system_class'],par_dict['maxeig'],par_dict['new_index'] = number_steadystates,steadystate_values_ss_n,ss_class,system_class,maxeig,[parID,ss_n]
                output_df = pd.concat([output_df,pd.DataFrame([par_dict], columns=par_dict.keys())], ignore_index=True)
        else:
            par_dict['ss_n'],par_dict['ss_list'],par_dict['ss_class'],par_dict['system_class'],par_dict['maxeig'],par_dict['new_index'] = 0, np.nan, np.nan,'no steady state', np.nan,[parID,0]
            output_df = pd.concat([output_df,pd.DataFrame([par_dict], columns=par_dict.keys())], ignore_index=True)

    output_df = output_df.set_index('new_index')
    return output_df
#Turing analysis carried out on a single parameter combination. The input is a dictionary with the corresponding parameters.
def detailed_turing_analysis_dict(par_dict, circuit_n,n_species,top_dispersion=5000,calculate_unstable=False):
    print(par_dict)
    steadystatelist, number_steadystates = findsteadystates(par_dict, circuit_n,n_species,n_initial_conditions=100) #input a dictionary with the parameters and returns (1) a list with the steady states and (2) the number of steady states.

    maxeig_list = []
    system_class_list = []
    ss_class_list = []
    eigenvalues_list=[]
    maxeig_list = []


    if number_steadystates > 0:
        for ss_n in range(number_steadystates): #perform linear stability analysis on all steady states found
            steadystate_values_ss_n = steadystatelist[ss_n]
            ss_class, system_class, eigenvalues, maxeig= dispersionrelation(par_dict,steadystate_values_ss_n, circuit_n,top_dispersion)

            system_class_list.append(system_class)
            ss_class_list.append(system_class)
            eigenvalues_list.append(eigenvalues)
            maxeig_list.append(maxeig)

    else:
        eigenvalues=[]
    return steadystatelist, number_steadystates, ss_class_list, system_class_list, eigenvalues_list, maxeig_list