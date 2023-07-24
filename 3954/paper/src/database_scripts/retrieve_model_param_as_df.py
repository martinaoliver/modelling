#%%
import sys
import os

pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
sys.path.append(modellingpath + '/lib')

from database.databaseFunctions import *
import pickle

#############

nsr_list=[0.01,0.05,0.1,0.2]
for nsr in nsr_list:




    #%%
    # Specify name of circuit and variant investigated
    circuit_n='14'
    parID=4187715
    # nsr=0.01
    old_variant='fitted7'
    variant=f'{old_variant}_gaussian{parID}_nsr{nsr}' #variant='fitted7'
    # Specifiy number of parameter sets in parameterset file to be loaded
    n_samples = 2000 #n_samples = 13700000


    model_param_dict = {'circuit_n': circuit_n, 'variant': variant, 'n_samples': n_samples}


    #%%
    result_df = query_modelParam_df_from_sql(model_param_dict)


    # pickle.dump(result_df, open('file.txt', 'wb' ) )
    pickle.dump(result_df, open(modellingpath + '/3954/paper/input/gaussian_parameterfiles/df_circuit%s_variant%s_%rparametersets.pkl'%(circuit_n,variant,n_samples), 'wb' ) )
    # with open(modellingpath + '/3954/paper/input/gaussian_parameterfiles/df_circuit%r_variant%s_%rparametersets.pkl'%(circuit_n,variant,n_samples), 'wb' ) as f:
        # pickle.dump(object, result_df)

    # %%



