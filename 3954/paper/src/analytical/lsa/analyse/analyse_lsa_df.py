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


def pieChart_lsa(valueCounts_dict,title,percentageCounts_dict,log=True):
    colors_dict={'simple stable':'grey','simple unstable':'grey','complex unstable':'grey','no steady state':'grey','hopf':'peachpuff','turing I hopf':'peachpuff','turing I':'coral','turing I oscillatory':'coral'}
    # colors=['grey','grey','grey','grey','peachpuff','peachpuff','peachpuff','coral','coral','coral']
    labels = []
    sizes = []
    colors= []
    percentages = []
    for x, y in valueCounts_dict.items():
        percentage = percentageCounts_dict[x]
        
        labels.append(f'{x} {np.round(percentage,4)} %')
        sizes.append(y)
        colors.append(colors_dict[x])
    if log==True:
        sizes = np.log(sizes)
    plt.pie(sizes,colors=colors, labels=labels )
    plt.tight_layout()

    plt.axis('equal')
    plt.title(title)


def hist_lsa(valueCounts_dict,title,percentageCounts_dict,log=True):
    colors_dict={'simple stable':'grey','simple unstable':'grey','complex unstable':'grey','no steady state':'grey','hopf':'peachpuff','turing I hopf':'peachpuff','turing I':'coral','turing I oscillatory':'coral'}
    # colors=['grey','grey','grey','grey','peachpuff','peachpuff','peachpuff','coral','coral','coral']
    labels = []
    sizes = []
    colors= []
    percentages = []
    for x, y in valueCounts_dict.items():
        percentage = percentageCounts_dict[x]
        
        labels.append(f'{x} {np.round(percentage,4)} %')
        sizes.append(y)
        colors.append(colors_dict[x])
    if log==True:
        sizes = np.log(sizes)
    plt.bar(sizes,colors=colors, labels=labels )
    plt.tight_layout()

    plt.axis('equal')
    plt.title(title)

    # plt.show()
#%%

# Specify name of circuit and variant investigated
circuit_n='circuit14'
variant='fitted7'
balance = 'balancedSemiBalanced'
# Specifiy number of parameter sets in parameterset file to be loaded
n_param_sets = 13700000#1000000
# Kce=100
print(f'{circuit_n}, Variant:{variant}, {balance}')
# lsa_df_circuit14_variantfitted7_13700000parametersets_balancedSemiBalanced.pkl
# df= pickle.load( open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/all_dataframes/lsa_df_%s_variant%s_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), "rb"))
# df= pickle.load( open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/all_dataframes/lsa_df_%s_variant%s_%rparametersets_balancedSemiBalanced.pkl'%(circuit_n,variant,n_param_sets), "rb"))
df= pickle.load( open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/all_dataframes/lsa_df_%s_variant%s_%rparametersets_%s.pkl'%(circuit_n,variant,n_param_sets,balance), "rb"))
# df= pickle.load( open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/all_dataframes/lsa_df_%s_variant%s_%rparametersets_%s_Kce%s.pkl'%(circuit_n,variant,n_param_sets,balance,Kce), "rb"))
# df= pickle.load( open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/all_dataframes/lsa_df_%s_variant%s_%rparametersets_balancedSemiBalanced.pkl'%(circuit_n,variant,n_param_sets), 'rb'))

print(df['system_class'].value_counts())

# %%

valueCounts_dict = dict(df['system_class'].value_counts())
percentageCounts_dict= {k: v/len(df)*100 for k, v in valueCounts_dict.items()}
title = f'{circuit_n} Variant {variant} {balance}, Kce {Kce}'
# plt.savefig(modellingpath + f'/3954/paper/out/analytical/pyPlots/piechart_{circuit_n}_variant{variant}_nsamples{n_param_sets}_{balance}_Kce{Kce}.pdf')
# plt.savefig(modellingpath + f'/3954/paper/out/analytical/pyPlots/piechart_{circuit_n}_variant{variant}_nsamples{n_param_sets}_{balance}_Kce{Kce}.png')
plt.bar(percentageCounts_dict.keys(), percentageCounts_dict.values()) 
plt.xticks(rotation=45, ha='right')
plt.tight_layout()

plt.yscale('log')
plt.show()
# dfunstable = df[df['system_class']=='simple unstable']




fig, ax = plt.subplots()
bars = ax.barh(list(percentageCounts_dict.keys()), np.round(list(percentageCounts_dict.values()),3))

ax.bar_label(bars)

for bars in ax.containers:
    ax.bar_label(bars)
ax.set_xscale('log')
plt.tight_layout()

plt.savefig(modellingpath + f'/3954/paper/out/analytical/pyPlots/barchart_{circuit_n}_variant{variant}_nsamples{n_param_sets}_{balance}_Kce{Kce}.pdf')
plt.savefig(modellingpath + f'/3954/paper/out/analytical/pyPlots/barchart_{circuit_n}_variant{variant}_nsamples{n_param_sets}_{balance}_Kce{Kce}.png')

plt.show()



#%%


valueCounts_dict = dict(df['system_class'].value_counts())
percentageCounts_dict= {k: v/len(df)*100 for k, v in valueCounts_dict.items()}
title = f'{circuit_n} Variant {variant} {balance}, Kce {Kce}'
pieChart_lsa(valueCounts_dict,title,percentageCounts_dict )
plt.savefig(modellingpath + f'/3954/paper/out/analytical/pyPlots/piechart_{circuit_n}_variant{variant}_nsamples{n_param_sets}_{balance}_Kce{Kce}.pdf')
plt.savefig(modellingpath + f'/3954/paper/out/analytical/pyPlots/piechart_{circuit_n}_variant{variant}_nsamples{n_param_sets}_{balance}_Kce{Kce}.png')
plt.show()
dfunstable = df[df['system_class']=='simple unstable']

# %%
#%%
from tqdm import tqdm
def weightTuring(dfLoc):
    for i in range(10):
        turing_lenght = len(dfLoc.loc[dfLoc['system_class'].isin(['turing I oscillatory','turing I','turing I hopf'])])
        # print(turing_lenght/ len(dfLoc))
        return turing_lenght/ len(dfLoc)

for i in tqdm(range(len(df))):
    # print(df.loc[i]['system_class'])
    dfLoc = df.loc[i]
    weight = weightTuring(dfLoc)
    df.at[i,'weightTuring']=weight
    # print('----')
    


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
    pickle.dump( instabilities_df, open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/instabilities_dataframes/instability_df_%s_variant%s_%rparametersets_%s_Kce%s.pkl'%(circuit_n,variant,n_param_sets,balance, Kce), "wb" ) )
    len(instabilities_df)
# instabilityComplex_df = pd.concat([complex_df, instabilities_df])
# pickle.dump( instabilityComplex_df, open(modellingpath + '/growth/out/analytical/instabilityComplex/instabilityComplex_df_%s_variant%s_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), "wb" ) )
# #values that have turing
saveTuring = True
if saveTuring == True:
    turingStates = ['turing I','turing I oscillatory', 'turing I hopf']  
    turing_df = df.loc[df['system_class'].isin(turingStates)]
    # turing_df = turing_df.xs(0, level=1)
    turing_df.index  = turing_df.index.droplevel(-1)
    # turing_df = turing_df.loc[turing_df['ss_n']==1]
    # turing_df = turing_df.sort_values(by=['maxeig'],  ascending=False)
    # pickle.dump( turing_df, open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/turing_dataframes/turing_df_%s_variant%s_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), "wb" ) )
    pickle.dump( turing_df, open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/turing_dataframes/turing_df_%s_variant%s_%rparametersets_%s_Kce%s.pkl'%(circuit_n,variant,n_param_sets,balance,Kce), "wb" ) )
#%%
compare_balance_dfs_systemclass=True
if compare_balance_dfs_systemclass == True:
    df0= pickle.load( open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/turing_dataframes/turing_df_%s_variant%s_%rparametersets_%s.pkl'%(circuit_n,variant,n_param_sets,'balanced'), "rb"))
    df1= pickle.load( open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/turing_dataframes/turing_df_%s_variant%s_%rparametersets_%s.pkl'%(circuit_n,variant,n_param_sets,'semiBalanced'), "rb"))
    df2= pickle.load( open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/turing_dataframes/turing_df_%s_variant%s_%rparametersets_%s.pkl'%(circuit_n,variant,n_param_sets,'notBalanced'), "rb"))


    graph_df = df0['system_class'].value_counts().rename('Original parameters (cir2var0)').to_frame()\
                .join(df1['system_class'].value_counts().rename('Converted parameters (cir2var2)').to_frame()).join(df2['system_class'].value_counts().rename('Converted parameters (cir2var3)').to_frame())
    graph_df.plot(kind='bar',figsize=(8, 4), color=['darkcyan','lightcoral', 'blue'], log=True)

    plt.legend()
    plt.title('Robustness of original params vs converted params')
    plt.ylabel('Number of Occurrences (log)', fontsize=12)
    plt.xlabel('System class', fontsize=12)
    plt.show()



# %%
#%%
compare_balance_dfs_general=False
if compare_balance_dfs_general == True:
    df0= pickle.load( open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/all_dataframes/lsa_df_%s_variant%s_%rparametersets_%s.pkl'%(circuit_n,variant,n_param_sets,'balanced'), "rb"))
    df1= pickle.load( open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/all_dataframes/lsa_df_%s_variant%s_%rparametersets_%s.pkl'%(circuit_n,variant,n_param_sets,'semiBalanced'), "rb"))
    df2= pickle.load( open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/all_dataframes/lsa_df_%s_variant%s_%rparametersets_%s.pkl'%(circuit_n,variant,n_param_sets,'notBalanced'), "rb"))

    

    df0_turing= pickle.load( open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/turing_dataframes/turing_df_%s_variant%s_%rparametersets_%s.pkl'%(circuit_n,variant,n_param_sets,'balanced'), "rb"))
    df1_turing= pickle.load( open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/turing_dataframes/turing_df_%s_variant%s_%rparametersets_%s.pkl'%(circuit_n,variant,n_param_sets,'semiBalanced'), "rb"))
    df2_turing= pickle.load( open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/turing_dataframes/turing_df_%s_variant%s_%rparametersets_%s.pkl'%(circuit_n,variant,n_param_sets,'notBalanced'), "rb"))

    robustness0 = len(df0_turing)/len(df0)*100
    robustness1 = len(df1_turing)/len(df1)*100
    robustness2 = len(df2_turing)/len(df2)*100

    plt.bar(['balanced', 'semiBalanced', 'notBalanced'],[robustness0, robustness1, robustness2])
    # graph_df = df0['system_class'].value_counts().rename('Original parameters (cir2var0)').to_frame()\
    #             .join(df1['system_class'].value_counts().rename('Converted parameters (cir2var2)').to_frame()).join(df2['system_class'].value_counts().rename('Converted parameters (cir2var3)').to_frame())
    # graph_df.plot(kind='bar',figsize=(8, 4), color=['darkcyan','lightcoral', 'blue'], log=True)

    plt.legend()
    plt.title('Robustness of circuit balancing')
    plt.ylabel('Percentage of Turing occurrences', fontsize=12)
    plt.tight_layout()
    plt.savefig(modellingpath + '/3954/paper/out/analytical/pyPlots/balancing_robustness_df_%s_variant%s_%rparametersets.png'%(circuit_n,variant,n_param_sets))
    plt.savefig(modellingpath + '/3954/paper/out/analytical/pyPlots/balancing_robustness_df_%s_variant%s_%rparametersets.pdf'%(circuit_n,variant,n_param_sets))
    plt.show()

# %%
