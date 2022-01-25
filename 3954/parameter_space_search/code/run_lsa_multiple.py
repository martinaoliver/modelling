##########################
#########IMPORTS##########
##########################

import os
import os.path
import sys
sixeqpath ='/rds/general/user/mo2016/home/Documents/modelling/3954'
modulepath = sixeqpath + '/modules'
sys.path.append(modulepath)
while True:
    try:
        from linear_stability_analysis import *
        break
    except ImportError:
        sixeqpath ='/Volumes/mo2016/home/Documents/modelling/3954/'
        modulepath = sixeqpath + '/modules'
        sys.path.append(modulepath)
        from linear_stability_analysis import *

import sys
import time
import pickle

#######################
#########CODE##########
#######################


circuit_n = 2 #ID of circuit we want to analyse
#(parameter sets provided correspond to circuit2 which is the one currently being implemented experimentally)
parID = 0 #takes the first parameter set of the dataframe... can choose any
n_species=6 #number of molecular species in circuit_n (#Circuit2 has 6 molecular species)
variant=0
n_param_sets=1000000
n_param_sets=10
start_batch_index = int(sys.argv[1])
batch_size = 3
#obtain a dictionary with some parameters to use in our analysis
df= pickle.load( open('../parameterfiles/df_circuit%r_variant%r_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), "rb" ) )
df = df.iloc[start_batch_index:start_batch_index+batch_size]
print('df loaded')
#Run analysis on 1M parameter sets
output_df = big_turing_analysis_df(df,circuit_n,n_species,print_parID=True)
print(output_df)
# pickle.dump(output_df, open('../results/output_dataframes/lsa_df_circuit%r_variant%r_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), 'wb'))
