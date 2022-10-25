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


circuit_n='circuit12'
variant=0

parID = 7 #takes the first parameter set of the dataframe... can choose any
n_species=9 #number of molecular species in circuit_n (#Circuit2 has 6 molecular species)

df_lenght = 1000
n_param_sets = 1000

start_batch_index = 0
batch_size = 1000
# obtain a dictionary with some parameters to use in our analysis
df= pickle.load( open(modellingpath + "/3954/paper/input/parameterfiles/df_%s_variant%r_%rparametersets.pkl"%(circuit_n,variant,n_param_sets), "rb"))
par_dict = df.loc[parID].to_dict()
#Run analysis on 1M parameter sets
# output_df = big_turing_analysis_df(df_batch,circuit_n,n_species,print_parID=False, tqdm_disable=False)
out = detailed_turing_analysis_dict(par_dict, circuit_n,n_species,top_dispersion=5000,calculate_unstable=False,steadystate=False)
plot_all_dispersion(out[4][2],9, crop=40)
print(out[3])
# multiple_df= pickle.load( open(modellingpath + "/growth/out/analytical/lsa_dataframes/lsa_df_%s_variant%r_%rparametersets.pkl"%(circuit_n,variant,n_param_sets), "rb"))
# par_dict = multiple_df.loc[parID].to_dict() #converts a dataframe row into a dictionary outputing a dictionary for a specific parameter set
# print(par_dict)