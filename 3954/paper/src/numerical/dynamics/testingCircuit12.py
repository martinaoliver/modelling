#############
###paths#####
#############
import sys
import os

pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
sys.path.append(modellingpath + '/lib')
#############
import pickle as pkl
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

circuit_n=12
variant= 0
nsamples=1000

df= pkl.load( open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/all_dataframes/lsa_df_circuit%s_variant%s_%rparametersets.pkl'%(circuit_n,variant,nsamples), "rb" ) )
print(df['ss_list'][0])
ss_list_variables = ['ssU','ssV','ssA','ssB','ssC','ssD','ssE','ssF', 'ssaTc']
df[ss_list_variables] = pd.DataFrame(df.ss_list.tolist(), index= df.index)
df
#plot distributions of steady states
for ss_variable in ss_list_variables:
    sns.displot(df[ss_variable],bins=5)
    plt.show()