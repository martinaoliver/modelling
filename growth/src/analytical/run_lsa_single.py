

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


circuit_n='turinghill'
variant=2 
n_species=2 #number of molecular species in circuit_n (#Circuit2 has 6 molecular species)

n_param_sets = 100000

# obtain a dictionary with some parameters to use in our analysis
df= pickle.load( open(modellingpath + "/growth/input/parameterfiles/df_%s_variant%r_%rparametersets.pkl"%(circuit_n,variant,n_param_sets), "rb"))
parID=32095
par_dict = df.loc[parID].to_dict() #converts a dataframe row into a dictionary outputing a dictionary for a specific parameter set
#Run analysis on 1M parameter sets
out = detailed_turing_analysis_dict(par_dict, circuit_n, n_species)
plot_all_dispersion(out[4][0],2, crop=100)
print(out[3])
multiple_df= pickle.load( open(modellingpath + "/growth/out/analytical/lsa_dataframes/lsa_df_%s_variant%r_%rparametersets.pkl"%(circuit_n,variant,n_param_sets), "rb"))
par_dict = multiple_df.loc[parID].to_dict() #converts a dataframe row into a dictionary outputing a dictionary for a specific parameter set
print(par_dict)