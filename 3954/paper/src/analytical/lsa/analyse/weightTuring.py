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
import matplotlib.pyplot as plt
import numpy as np
import pickle as pkl

def pieChart_lsa(valueCounts_dict,title,log=True):
    colors_dict={'simple stable':'grey','simple unstable':'grey','complex unstable':'grey','no steady state':'grey','hopf':'peachpuff','turing I hopf':'peachpuff','turing I':'coral','turing I oscillatory':'coral'}
    # colors=['grey','grey','grey','grey','peachpuff','peachpuff','peachpuff','coral','coral','coral']
    labels = []
    sizes = []
    colors= []
    for x, y in valueCounts_dict.items():
        labels.append(x)
        sizes.append(y)
        colors.append(colors_dict[x])
    if log==True:
        sizes = np.log(sizes)
    plt.pie(sizes,colors=colors, labels=labels)
    plt.axis('equal')
    plt.title(title)
    plt.show()
#%%

# Specify name of circuit and variant investigated
circuit_n='circuit14'
variant='2nd'
# Specifiy number of parameter sets in parameterset file to be loaded
n_param_sets = 1000000
balance='notBalanced'
print(f'Circuit:{circuit_n}, Variant:{variant}')

df= pickle.load( open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/all_dataframes/lsa_df_%s_variant%s_%rparametersets_%s.pkl'%(circuit_n,variant,n_param_sets,balance), "rb"))
print(df['system_class'].value_counts())


#%%
from tqdm import tqdm
def weightTuring(dfLoc):
    turing_lenght = len(dfLoc.loc[dfLoc['system_class'].isin(['turing I oscillatory','turing I', 'turing I hopf'])])
        # turing_
        # lenght = len(dfLoc.loc[dfLoc['system_class'].isin(['simple stable','turing I'])])
        # print(turing_lenght/ len(dfLoc))
    return turing_lenght/ len(dfLoc)

for i in tqdm(range(n_param_sets)):
    # print(df.loc[i]['system_class'])
    dfLoc = df.loc[i]
    weight = weightTuring(dfLoc)
    df.at[i,'weightTuring']=weight
print('Turing robustness',np.sum(df['weightTuring'])/n_param_sets)

pkl.dump(df, open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/all_dataframes/lsa_df_%s_variant%s_%rparametersets_%s_weightTuring.pkl'%(circuit_n,variant,n_param_sets,balance), "wb"))
# %%
dfBalanced= pickle.load( open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/all_dataframes/lsa_df_%s_variant%s_%rparametersets_balanced_weightTuring.pkl'%(circuit_n,variant,n_param_sets), "rb"))
semiBalanced= pickle.load( open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/all_dataframes/lsa_df_%s_variant%s_%rparametersets_semiBalanced_weightTuring.pkl'%(circuit_n,variant,n_param_sets), "rb"))
notBalanced= pickle.load( open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/all_dataframes/lsa_df_%s_variant%s_%rparametersets_notBalanced_weightTuring.pkl'%(circuit_n,variant,n_param_sets), "rb"))


# %%
import seaborn as sns
nsamples=1000000
balancedRobustness = np.sum(dfBalanced['weightTuring'])/nsamples*100
semiBalancedRobustness = np.sum(semiBalanced['weightTuring'])/nsamples*100
notBalancedRobustness = np.sum(notBalanced['weightTuring'])/nsamples*100
plt.bar(['Balanced', 'semi balanced', 'not balanced'],[balancedRobustness,semiBalancedRobustness,notBalancedRobustness], color=['green','yellow', 'red'])
plt.ylim(0,1e-2)
plt.ylabel('Turing robustness %')



# %%
