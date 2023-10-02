
#%%
###paths#####
#############
import sys
import os

from importlib_metadata import distribution
pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
sys.path.append(modellingpath + '/lib')
#############

from analytical.linear_stability_analysis import big_turing_analysis_df, detailed_turing_analysis_dict
from randomfunctions import plot_all_dispersion

import pickle

#######################
#########CODE##########
#######################

circuit_n='turinghill'; n_species=2

df= pickle.load( open(modellingpath + f'/growth/out/analytical/instability/instability_df_circuitturinghill_variant8-9_combinedparametersets.pkl', 'rb'))
df.set_index('ssID', append=True, inplace=True)
#%%
parID=640178; ssID=0
par_dict = df.loc[(parID,ssID)].to_dict()

out = detailed_turing_analysis_dict(par_dict, circuit_n, n_species)
steadystatelist, number_steadystates, ss_class_list, system_class_list, eigenvalues_list, maxeig_list, complex_dispersion_list = out
plot_all_dispersion(out[4][0],2, crop=100)
print(out[3])


