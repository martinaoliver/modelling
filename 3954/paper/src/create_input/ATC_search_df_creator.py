#############
###paths#####
#############
import sys
import os




pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1]
sys.path.append(modellingpath + '/lib')
#############
from equations.parameterCreation_functions import *
#############
import numpy as np
import pandas as pd
import pickle 
import matplotlib.pyplot as plt


#%%

# Specify name of circuit and variant investigated
circuit_n=14
variant='2nd'
balance = 'balanced'
# Specifiy number of parameter sets in parameterset file to be loaded
n_samples = 1000000

print(f'Circuit:{circuit_n}, Variant:{variant}')
df= pickle.load( open(modellingpath + "/3954/paper/input/balanced_parameterfiles/df_circuit%s_variant%s_%rparametersets_%s.pkl"%(circuit_n,variant,n_samples, balance), "rb"))

kce_values_list = [0.1,1,10,100,1000]

for Kce in kce_values_list:
    df['Kce'] = Kce
    plt.hist(df['Kce'])
    plt.show()
    pickle.dump(df, open(modellingpath + '/3954/paper/input/balanced_parameterfiles/df_circuit%s_variant%s_%rparametersets_balanced_kce%s.pkl'%(circuit_n,variant,n_samples,Kce), 'wb'))
