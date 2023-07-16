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

#%%
#############

# #%%
# Specify name of circuit and variant investigated

circuit_n='turinghill'
variant=9

# Specifiy number of parameter sets in parameterset file to be loaded
n_samples = 2000000
# n_samples = 10

print(f'Circuit:{circuit_n}, Variant:{variant}')
lsa_df = pickle.load( open(modellingpath + '/growth/out/analytical/lsa_dataframes/lsa_df_%s_variant%s_%rparametersets.pkl'%(circuit_n,variant,n_samples), "rb"))
lsa_df_cropped = lsa_df.iloc[:10]





#%%
analyticalOutput_df_to_sql(lsa_df_cropped, circuit_n, variant, n_samples)

# %%
