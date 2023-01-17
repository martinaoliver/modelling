#############

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


circuit_n='circuit14'
variant='0nd'

parID = 7 #takes the first parameter set of the dataframe... can choose any
n_species=6 #number of molecular species in circuit_n (#Circuit2 has 6 molecular species)

df_lenght = 10
n_param_sets = 10

# obtain a dictionary with some parameters to use in our analysis
df= pickle.load( open(modellingpath + "/3954/paper/input/lhs_parameterfiles/df_%s_variant%s_%rparametersets.pkl"%(circuit_n,variant,n_param_sets), "rb"))
par_dict = df.loc[parID].to_dict()
#Run analysis on 1M parameter sets
# output_df = big_turing_analysis_df(df_batch,circuit_n,n_species,print_parID=False, tqdm_disable=False)
out = detailed_turing_analysis_dict(par_dict, circuit_n,n_species,top_dispersion=5000,calculate_unstable=False,steadystate=False)
plot_all_dispersion(out[4][2],n_species, crop=30)

print(out[3])

