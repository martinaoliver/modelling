#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 12:33:13 2020

@author: mo2016
"""


from findsteadystates_functions import findsteadystates
from dispersionrelation_functions import dispersionrelation


#Turing analysis carried out on a dataframe. The input is a df with every parameter set.
def big_turing_analysis_df(df,circuit_n,n_species):
    len_df = len(df) #lenght of dataframe (number of parameter sets to analyse)

    for n in range(len_df):

        print ('')
        print ('Parameter set %r ' %n)
        par_dict = df.iloc[n].to_dict() #converts a dataframe row into a dictionary outputing a dictionary for a specific parameter set
        steadystatelist, number_steadystates = findsteadystates(par_dict,circuit_n,n_species) #input a dictionary with the parameters and returns (1) a list with the steady states and (2) the number of steady states.

        if number_steadystates > 0:
            for ss_n in range(number_steadystates): #perform linear stability analysis on all steady states found
                steadystate_values_ss_n = steadystatelist[ss_n]
                print (steadystate_values_ss_n)
                maxeig, pattern_class, eigenvalues = dispersionrelation(par_dict,steadystate_values_ss_n, circuit_n)
        else:
            print('no steady states')

#Turing analysis carried out on a single parameter combination. The input is a dictionary with the corresponding parameters.
def detailed_turing_analysis_dict(par_dict, circuit_n,n_species):

    steadystatelist, number_steadystates = findsteadystates(par_dict, circuit_n,n_species) #input a dictionary with the parameters and returns (1) a list with the steady states and (2) the number of steady states.

    maxeig_list = []
    pattern_class_list = []
    eigenvalues_list=[]

    if number_steadystates > 0:
        for ss_n in range(number_steadystates): #perform linear stability analysis on all steady states found

            steadystate_values_ss_n = steadystatelist[ss_n]
            print (steadystate_values_ss_n)
            maxeig, pattern_class, eigenvalues = dispersionrelation(par_dict,steadystate_values_ss_n, circuit_n)
            maxeig_list.append(maxeig)
            pattern_class_list.append(pattern_class)
            eigenvalues_list.append(eigenvalues)

    else:
        eigenvalues=[]
    return steadystatelist, number_steadystates, maxeig_list, pattern_class_list, eigenvalues_list
