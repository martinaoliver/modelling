

#%%
#############
###paths#####
#############
import sys
import os

pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
sys.path.append(modellingpath + '/lib')
#############


import pickle
import pandas as pd


from numpy import random
random.seed(1)
#execution parameters
circuit_n=14
variant='fitted7'
balance = 'balancedSemiBalanced'

n_samples = 13700000
n_param_sets_new = 2000
# original_df= pickle.load( open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/all_dataframes/lsa_df_circuit%s_variant%s_%rparametersets_%s.pkl'%(circuit_n,variant,n_samples,balance), "rb"))
original_df= pickle.load( open(modellingpath + '/3954/paper/input/fitted_parameterfiles/df_circuit%s_variant%s_%rparametersets_%s.pkl'%(circuit_n,variant,n_samples,balance), "rb"))
#drop balanced column
original_df  = original_df.drop(['balance'], axis=1)
parID = 4187715
# nsr_list=[0.01,0.02,0.04, 0.06,0.08, 0.1,0.21]
nsr_list=[0.01,0.05,0.1,0.2]
# nsr_list=[0.05,0.1]





#%%
for nsr in nsr_list:
    # open parameter dictionaries
    series = original_df.loc[parID]
    
    df = pd.DataFrame()
    series = df.append(series)
    # df = pd.concat([series]*30000, ignore_index=True)
    df = pd.concat([series]*n_param_sets_new, ignore_index=True)
    

    #ADD NORMAL DISTRIBUTIONS
    n_list = ['nvd', 'nub', 'nda', 'nfe','nee', 'neb', 'nce']
    parameter_list = [parameter for parameter in df.columns]
    for parameter in parameter_list:
        if parameter not in n_list:
            df[parameter] = random.normal(loc=df[parameter].iloc[0], scale=df[parameter].iloc[0]*nsr, size=len(df))
    # n_list = ['nvd', 'nub', 'nda', 'nfe','nee', 'neb', 'nce']

    # for n in n_list:
    #     df[n] = original_df[n]

    # print((df < 0).values.any())
    new_variant=f'{variant}_gaussian{parID}_nsr{nsr}'
    pickle.dump( df, open(modellingpath + '/3954/paper/input/gaussian_parameterfiles/df_circuit%r_variant%s_%rparametersets.pkl'%(circuit_n,new_variant,n_param_sets_new), "wb" ) )
    # print(nsr)

print('df',df)

# %%

