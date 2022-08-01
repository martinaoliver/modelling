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


import pickle

#######################
#########CODE##########
#######################


circuit_n='schnakenberg'
variant=0 

parID = 0 #takes the first parameter set of the dataframe... can choose any
n_species=2 #number of molecular species in circuit_n (#Circuit2 has 6 molecular species)

n_param_sets=10

start_batch_index = 0
batch_size = 100
#obtain a dictionary with some parameters to use in our analysis
df= pickle.load( open(modellingpath + "/growth/input/parameterfiles/df_%s_variant%r_%rparametersets.pkl"%(circuit_n,variant,n_param_sets), "rb"))
print(len(df))
df = df.iloc[start_batch_index:start_batch_index+batch_size]
print('df loaded')
#Run analysis on 1M parameter sets
output_df = big_turing_analysis_df(df,circuit_n,n_species,print_parID=True)
print(output_df)
print(output_df.columns)
pickle.dump(output_df, open(modellingpath + '/growth/out/lsa_dataframes/lsa_df_%s_variant%r_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), 'wb'))