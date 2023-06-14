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
variant='fitted6'
# balance = 'balanced'
# Specifiy number of parameter sets in parameterset file to be loaded
n_param_sets = 8000000

print(f'Circuit:{circuit_n}, Variant:{variant}')

# df= pickle.load( open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/all_dataframes/lsa_df_%s_variant%s_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), "rb"))
# df= pickle.load( open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/all_dataframes/lsa_df_%s_variant%s_%rparametersets_balancedSemiBalanced.pkl'%(circuit_n,variant,n_param_sets), "rb"))
# df= pickle.load( open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/all_dataframes/lsa_df_%s_variant%s_%rparametersets_%s.pkl'%(circuit_n,variant,n_param_sets,balance), "rb"))
df= pickle.load( open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/all_dataframes/lsa_df_%s_variant%s_%rparametersets_balancedSemiBalanced.pkl'%(circuit_n,variant,n_param_sets), 'rb'))

print(df['system_class'].value_counts())


#%%
from tqdm import tqdm
def weightTuring(dfLoc):
    for i in range(10):
        turing_lenght = len(dfLoc.loc[dfLoc['system_class'].isin(['turing I oscillatory','turing I'])])
        # print(turing_lenght/ len(dfLoc))
        return turing_lenght/ len(dfLoc)

for i in tqdm(range(1000000)):
    # print(df.loc[i]['system_class'])
    dfLoc = df.loc[i]
    weight = weightTuring(dfLoc)
    df.at[i,'weightTuring']=weight
    # print('----')
    


# %%

valueCounts_dict = dict(df['system_class'].value_counts())
title = f'{circuit_n} Variant {variant} Not-Balanced'
pieChart_lsa(valueCounts_dict,title)
dfunstable = df[df['system_class']=='simple unstable']
#%%
if balanced == True:
    

# #values for which complex dispersion = true
# complex_df = df[df['complex_dispersion']==True]
# complex_df.index  = complex_df.index.droplevel(-1)
#%%
# #values that have instabilities
saveInstabilities = True
if saveInstabilities ==True:
    instabilities = ['turing I', 'turing II', 'turing I hopf', 'turing I oscillatory', 'turing II hopf','hopf', 'turing semi-hopf']  
    instabilities_df = df.loc[df['system_class'].isin(instabilities)]
    instabilities_df.index  = instabilities_df.index.droplevel(-1)
    print(instabilities_df['system_class'].value_counts())
    pickle.dump( instabilities_df, open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/instabilities_dataframes/instability_df_%s_variant%s_%rparametersets_%s.pkl'%(circuit_n,variant,n_param_sets,balance), "wb" ) )
    len(instabilities_df)
# instabilityComplex_df = pd.concat([complex_df, instabilities_df])
# pickle.dump( instabilityComplex_df, open(modellingpath + '/growth/out/analytical/instabilityComplex/instabilityComplex_df_%s_variant%s_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), "wb" ) )

# #values that have turing
saveTuring = True
if saveTuring == True:
    turingStates = ['turing I','turing I oscillatory']  
    turing_df = df.loc[df['system_class'].isin(turingStates)]
    # turing_df = turing_df.xs(0, level=1)
    turing_df.index  = turing_df.index.droplevel(-1)
    # turing_df = turing_df.loc[turing_df['ss_n']==1]
    # turing_df = turing_df.sort_values(by=['maxeig'],  ascending=False)
    # pickle.dump( turing_df, open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/turing_dataframes/turing_df_%s_variant%s_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), "wb" ) )
    pickle.dump( turing_df, open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/turing_dataframes/turing_df_%s_variant%s_%rparametersets_%s.pkl'%(circuit_n,variant,n_param_sets,balance), "wb" ) )
#%%
compare_two_dfs=False
if compare_two_dfs == True:
    df0= pickle.load( open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/all_dataframes/lsa_df_%s_variant%s_%rparametersets.pkl'%(circuit_n,0,n_param_sets), "rb"))
    df2= pickle.load( open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/all_dataframes/lsa_df_%s_variant%s_%rparametersets.pkl'%(circuit_n,2,n_param_sets), "rb"))


    graph_df = df0['system_class'].value_counts().rename('Original parameters (cir2var0)').to_frame()\
                .join(df2['system_class'].value_counts().rename('Converted parameters (cir2var2)').to_frame())
    graph_df.plot(kind='bar',figsize=(8, 4), color=['darkcyan','lightcoral'], log=True)

    plt.legend()
    plt.title('Robustness of original params vs converted params')
    plt.ylabel('Number of Occurrences (log)', fontsize=12)
    plt.xlabel('System class', fontsize=12)
    plt.show()

# %%
