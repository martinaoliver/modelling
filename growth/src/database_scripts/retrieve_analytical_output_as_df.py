#%%
import sys
import os

pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
sys.path.append(modellingpath + '/lib')

from database.databaseFunctions import *
import pickle




circuit_n='\'turinghill\'';variant='\'9\'';n_species=6;n_samples = 2000000



model_param_dict = {'circuit_n': circuit_n, 'variant': variant, 'n_samples': n_samples}









#%%





result_df = query_analyticalOutput_df_from_sql(model_param_dict)
pickle.dump(result_df, open( modellingpath + f'/growth/out/analytical/lsa_dataframes/lsa_df_circuitturinghill_variant9_{n_samples}parametersets.pkl', "wb" ) )

print(result_df.columns)
# %%
