#%%

import sys
import os

pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
sys.path.append(modellingpath + '/lib')

from database.databaseFunctions import *
import pickle
import psycopg2
import time
import numpy as np
from tqdm import tqdm

#%%
#############

# #%%
# Specify name of circuit and variant investigated



# circuit_n='turinghill'
# variant=int(sys.argv[1])
# seed=int(sys.argv[2])
# n_species=2
# # Specifiy number of parameter sets in parameterset file to be loaded
# n_samples = int(sys.argv[3])

circuit_n='turinghill'
variant=9
seed=0
n_species=2
# Specifiy number of parameter sets in parameterset file to be loaded
n_samples = 1000000




lhs_df = pickle.load( open(modellingpath + '/growth/input/parameterfiles/df_%s_variant%s_%rparametersets_seed%s.pkl'%(circuit_n,variant,n_samples,seed), "rb"))
lsa_df = pickle.load( open(modellingpath + '/growth/out/analytical/lsa_dataframes/lsa_df_%s_variant%s_%rparametersets_seed%s.pkl'%(circuit_n,variant,n_samples, seed), "rb"))


#values that have instabilities

#%%
instabilities = ['turing I', 'turing II', 'turing I hopf', 'turing I oscillatory', 'turing II hopf','hopf', 'turing semi-hopf']  
# instabilities = ['complex unstable']
# instabilities = ['simple stable']  
lsa_df_instabilities = lsa_df.loc[lsa_df['system_class'].isin(instabilities)]
print('seed',seed)
print(lsa_df_instabilities)
lhs_df_instabilities = lhs_df.loc[lsa_df_instabilities.index.get_level_values(0)]
lhs_df_instabilities = lhs_df_instabilities.drop_duplicates()

# lhs_df_instabilities = lhs_df_instabilities.loc[:2000]
# lsa_df_instabilities = lsa_df_instabilities.loc[:2000]
#%%
modelParam_df_to_sql(lhs_df_instabilities, circuit_n, variant, n_samples)

analyticalOutput_df_to_sql(lsa_df_instabilities, circuit_n, variant, n_samples)






# %%
