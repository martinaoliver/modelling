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

circuit_n='turinghill'
variant=2

# Specifiy number of parameter sets in parameterset file to be loaded
n_samples = 10

print(f'Circuit:{circuit_n}, Variant:{variant}')
lhs_df = pickle.load( open(modellingpath + '/growth/input/parameterfiles/df_%s_variant%s_%rparametersets.pkl'%(circuit_n,variant,n_samples), "rb"))
lhs_df = lhs_df.iloc[:10]
modelParam_df_to_sql(lhs_df, circuit_n, variant, n_samples)
