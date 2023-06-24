#############

###paths#####
#############
import sys
import os

pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
sys.path.append(modellingpath + '/lib')
#############

from analytical.linear_stability_analysis import big_turing_analysis_df, detailed_turing_analysis_dict

from randomfunctions import *

import pickle

#######################
#########CODE##########
#######################








circuit_n='circuit14'
variant='fitted7'
# variant='2nd'
balance = 'balancedSemiBalanced'
parID = 6160552 #takes the first parameter set of the dataframe... can choose any
n_species=6 #number of molecular species in circuit_n (#Circuit2 has 6 molecular species)

n_samples = 13700000
# obtain a dictionary with some parameters to use in our analysis
# df= pickle.load( open(modellingpath + "/3954/paper/input/balanced_parameterfiles/df_%s_variant%s_%rparametersets_%s.pkl"%(circuit_n,variant,n_samples, balance), "rb"))
df= pickle.load( open(modellingpath + "/3954/paper/input/fitted_parameterfiles/df_%s_variant%s_%rparametersets.pkl"%(circuit_n,variant,n_samples), "rb"))
# df= pickle.load( open(modellingpath + "/3954/paper/input/fitted_parameterfiles/df_%s_variant%s_%rparametersets.pkl"%(circuit_n,variant,n_param_sets), "rb"))
# n_analysed_samples = 10
df = df.iloc[parID-1:parID+1]
par_dict = df.loc[parID].to_dict()
#Run analysis on 1M parameter sets
# output_df = big_turing_analysis_df(df,circuit_n,variant,n_samples, n_species,print_parID=False, tqdm_disable=False)
out = detailed_turing_analysis_dict(par_dict, circuit_n,n_species,top_dispersion=100,calculate_unstable=False,steadystate=False)
plot_all_dispersion(out[4][0],n_species, crop=100, top=100)
plt.show()
plot_highest_dispersion(out[4][0], crop=10000, top=100)
plt.title(f'Turing instability: {circuit_n}, {variant}, parID {parID}')
plt.show()

# print(out[3])


