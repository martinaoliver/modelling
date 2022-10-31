#############
###paths#####
#############
import sys
import os

from importlib_metadata import distribution
pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
sys.path.append(modellingpath + '/lib')
#############


import pickle
import pandas as pd


from numpy import random
random.seed(1)
#execution parameters
circuit_n=2
variant=0
nsamples = 1000000
nsamples_new = 2000
original_df= pickle.load( open(modellingpath + '/3954/paper/input/lhs_parameterfiles/df_circuit%s_variant%s_%rparametersets.pkl'%(circuit_n,variant,nsamples), "rb" ) )
parID = 48257
nsr_list=[0.01,0.02,0.04, 0.06,0.08, 0.1,0.21]
for nsr in nsr_list:
    # open parameter dictionaries
    series = original_df.loc[parID]
    df = pd.DataFrame()
    series = df.append(series)
    # df = pd.concat([series]*30000, ignore_index=True)
    df = pd.concat([series]*nsamples_new, ignore_index=True)

    #ADD NORMAL DISTRIBUTIONS
    parameter_list = [parameter for parameter in df.columns]
    for parameter in parameter_list:
        df[parameter] = random.normal(loc=df[parameter].iloc[0], scale=df[parameter].iloc[0]*nsr, size=len(df))
    n_list = [ 'nbd', 'nab', 'nda', 'nfe', 'nee', 'neb', 'nce']
    for n in n_list:
        df[n] = original_df[n]
    print((df < 0).values.any())
    pickle.dump( df, open(modellingpath + '/3954/paper/input/gaussian_parameterfiles/df_circuit%r_variant%sgaussian%snsr_%rparametersets.pkl'%(circuit_n,48257,nsr,nsamples_new), "wb" ) )
    print(nsr)
    print(df)

