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
variant= 4
n_species=2
# Specifiy number of parameter sets in parameterset file to be loaded
n_param_sets = 2000000


df= pickle.load( open(modellingpath + '/growth/out/analytical/turing/turing_df_%s_variant%r_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), "rb"))

sns.histplot(df['estimated_wvl'], bins=100, kde=True)
plt.xlabel('Estimated wavelenght (mm) from LSA')