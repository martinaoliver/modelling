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
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

print('heehehe')


circuit_n='turinghill'
variant= 11
n_species=2
# Specifiy number of parameter sets in parameterset file to be loaded
n_param_sets = 2000000


# df= pickle.load( open(modellingpath + '/growth/out/analytical/turing/turing_df_%s_variant%r_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), "rb"))
df= pickle.load( open(modellingpath + '/growth/out/analytical/lsa_dataframes/lsa_df_%s_variant%r_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), "rb"))
instabilities = ['turing I', 'turing II', 'turing I hopf', 'turing I oscillatory', 'turing II hopf', 'turing semi-hopf']  
instabilities_df = df.loc[df['system_class'].isin(instabilities)]
sns.histplot(df['estimated_wvl'], bins=100, kde=True)
plt.xlabel('Estimated wavelenght (mm) from LSA')
plt.show()