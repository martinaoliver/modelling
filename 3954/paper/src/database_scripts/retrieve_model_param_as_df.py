#%%
import sys
import os

pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
sys.path.append(modellingpath + '/lib')

from database.databaseFunctions import *
import pickle
#############


# #%%
# Specify name of circuit and variant investigated
circuit_n='circuit14'
variant='fitted7'
# Specifiy number of parameter sets in parameterset file to be loaded
n_samples = 13700000

print(f'Circuit:{circuit_n}, Variant:{variant}')
# lhs_df = pickle.load( open(modellingpath + '/3954/paper/input/lhs_parameterfiles/df_%s_variant%s_%rparametersets.pkl'%(circuit_n,variant,n_samples), "rb"))
lhs_df = pickle.load( open(modellingpath + '/3954/paper/input/fitted_parameterfiles/df_%s_variant%s_%rparametersets_balancedSemiBalanced.pkl'%(circuit_n,variant,n_samples), "rb"))
# lhs_df = pickle.load( open(modellingpath + '/3954/paper/input/fitted_parameterfiles/df_%s_variant%s_%rparametersets.pkl'%(circuit_n,variant,n_samples), "rb"))
batch_size = int(len(lhs_df/100))+1
lhs_df = lhs_df.iloc[:100000]
print(lhs_df.columns)
print('sgxjfh')

#%%
modelParam_df_to_sql(lhs_df, circuit_n, variant, n_samples)

# %%
