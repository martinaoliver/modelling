import sys
import os
pwd = os.getcwd()
root = pwd.split("home", 1)[0]
modelling_home = root + 'home/Documents/modelling'
modelling_ephemeral = root + 'ephemeral/Documents/modelling'
modulepath = modelling_home + '/3954/modules'
sys.path.append(modulepath)
import pickle
import pandas as pd
from numpy import random
random.seed(1)
#execution parameters
circuit=2
variant=0
n_parametersets = 10000
original_df = pd.read_pickle("df_circuit%r_variant%r_%rparametersets.pkl"%(circuit,variant,n_parametersets))
var_list=[0.01,0.02,0.04, 0.06,0.08, 0.1,0.23]
for var in var_list:
    # open parameter dictionaries
    series = original_df.loc[5716]
    df = pd.DataFrame()
    series = df.append(series)
    # df = pd.concat([series]*30000, ignore_index=True)
    df = pd.concat([series]*1000, ignore_index=True)

    #ADD NORMAL DISTRIBUTIONS
    parameter_list = ['Va','Vb','Vc','Vd','Ve','Vf','ba','bb','bc','bd','be','bf','d_B','kaa','kbd','kce','kda','keb','kee','kfe','mua','mulva']
    for parameter in parameter_list:
        df[parameter] = random.normal(loc=df[parameter].iloc[0], scale=df[parameter].iloc[0]*var, size=len(df))
    print((df < 0).values.any())
    pickle.dump( df, open( "5716gaussian/df_circuit%r_variant%s_%rparametersets_%rvar.pkl"%(circuit,'5716gaussian',1000,var), "wb" ) )
    print(var)
    print(df)