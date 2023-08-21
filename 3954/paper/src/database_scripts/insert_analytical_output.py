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
import math





#%%
#############

# #%%
# Specify name of circuit and variant investigated
circuit_n='14' #circuit_n='circuit14'

variant='fitted7'#'2nd' #'fitted7'#
balance='balancedSemiBalanced'
Kce=100
# Specifiy number of parameter sets in parameterset file to be loaded
n_samples = 13700000#1000000 #n_samples = 13700000#

print(f'Circuit:{circuit_n}, Variant:{variant}')
# lsa_df = pickle.load( open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/all_dataframes/lsa_df_circuit%s_variant%s_%rparametersets_%s_Kce%s.pkl'%(circuit_n,variant,n_samples,balance,Kce), "rb"))
lsa_df = pickle.load( open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/all_dataframes/lsa_df_circuit%s_variant%s_%rparametersets_%s.pkl'%(circuit_n,variant,n_samples,balance), "rb"))
#%%
# lsa_df_cropped = lsa_df.loc[12837401, 13699996]
# analyticalOutput_df_to_sql(lsa_df_cropped, circuit_n, variant, n_samples)
# lsa_df_cropped = lsa_df.iloc[:10]
# lsa_df = pickle.load(  open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/turing_dataframes/turing_df_circuit%s_variant%s_%sparametersets_balanced_Kce%s.pkl'%(circuit_n,variant,n_samples,Kce), "rb" ))


# lsa_df_cropped = lsa_df[lsa_df['ss_n']==0]
#%% 
#Insert in batches
batch_insert=True
if batch_insert==True:
    batch_size=10000
    batch_indices = list(range(0, len(lsa_df), batch_size))

    for n in tqdm(range(int(len(lsa_df)/batch_size))):
        print(n)
        lsa_df_cropped = lsa_df.iloc[batch_indices[n]:batch_indices[n] + batch_size+1 ]
        analyticalOutput_df_to_sql(lsa_df_cropped, circuit_n, variant, n_samples)
# #%%
# # lsa_df_cropped = lsa_df.iloc[:10]
# analyticalOutput_df_to_sql(lsa_df, circuit_n, variant, n_samples)
# # 
# %%
