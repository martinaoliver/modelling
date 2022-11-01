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

#############
###circuit2###
#############
print('circuit12')
circuit_n=12
variant= 0
nsamples=1000000

df= pkl.load( open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/all_dataframes/lsa_df_circuit%s_variant%s_%rparametersets.pkl'%(circuit_n,variant,nsamples), "rb" ) )
df = df.loc[:10000]
ss_list_variables = ['ssU','ssV','ssA','ssB','ssC','ssD','ssE','ssF', 'ssaTc']
df[ss_list_variables] = pd.DataFrame(df.ss_list.tolist(), index= df.index)
print(df.describe())

#plot distributions of steady states
fig,axs = plt.subplots(3,3,figsize=(15,15))
#flatten axes
axs = axs.flatten()
for count, ss_variable in enumerate(ss_list_variables):
    # df = df[df[ss_variable]< np.percentile(df[ss_variable],99.9)]
    # df = df[df[ss_variable]> np.percentile(df[ss_variable],0.1)]
    sns.histplot(df[ss_variable],bins=1000, log_scale=True, ax = axs[count])
plt.show()


#############
###circuit12###
#############
print('circuit2')
circuit_n=2
variant= 0
nsamples=1000000

df= pkl.load( open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/all_dataframes/lsa_df_circuit%s_variant%s_%rparametersets.pkl'%(circuit_n,variant,nsamples), "rb" ) )
df = df.loc[:10000]
# ss_list_variables = ['ssU','ssV','ssA','ssB','ssC','ssD','ssE','ssF', 'ssaTc']
ss_list_variables = ['ssA','ssB','ssC','ssD','ssE','ssF']
df[ss_list_variables] = pd.DataFrame(df.ss_list.tolist(), index= df.index)
# print(df['ssU'])

#plot distributions of steady states
for ss_variable in ss_list_variables:
    sns.displot(df[ss_variable],bins=1000, log_scale=True)
    plt.show()


print(df.describe())
