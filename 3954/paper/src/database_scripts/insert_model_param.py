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
# lhs_df = pickle.load( open(modellingpath + '/3954/paper/input/fitted_parameterfiles/df_%s_variant%s_%rparametersets_balancedSemiBalanced.pkl'%(circuit_n,variant,n_samples), "rb"))
lhs_df = pickle.load( open(modellingpath + '/3954/paper/input/fitted_parameterfiles/df_%s_variant%s_%rparametersets.pkl'%(circuit_n,variant,n_samples), "rb"))


batch_size=100000
batch_indices = list(range(0, len(lhs_df), batch_size))


#%%
for n in range (int(len(lhs_df)/batch_size)):
    print(n)
    lhs_df_cropped = lhs_df.iloc[batch_indices[n]:batch_indices[n] + batch_size+1 ]
    modelParam_df_to_sql(lhs_df_cropped, circuit_n, variant, n_samples)



# %%
