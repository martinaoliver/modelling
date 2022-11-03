###paths#####
import sys
import os
pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
sys.path.append(modellingpath + '/lib')
#############

from analytical.linear_stability_analysis import *


import pickle

#######################
#########CODE##########
#######################


circuit_n = 'circuit2' #ID of circuit we want to analyse
#(parameter sets provided correspond to circuit2 which is the one currently being implemented experimentally)
parID = 1 #takes the first parameter set of the dataframe... can choose any
n_species=6 #number of molecular species in circuit_n (#Circuit2 has 6 molecular species)
variant=1
n_param_sets=1000000
# n_param_sets=10
# start_batch_index = int(sys.argv[1])
start_batch_index = 0
batch_size = 3
#obtain a dictionary with some parameters to use in our analysis
df= pickle.load( open(modellingpath + '/3954/parameter_space_search/lhs_parameterfiles/df_%s_variant%r_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), "rb" ) )
print(len(df))
df = df.iloc[start_batch_index:start_batch_index+batch_size]
print('df loaded')
#Run analysis on 1M parameter sets
output_df = big_turing_analysis_df(df,circuit_n,n_species,print_parID=True)
print(output_df)
print(output_df.columns)
# pickle.dump(output_df, open('../results/output_dataframes/lsa_df_circuit%r_variant%r_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), 'wb'))
