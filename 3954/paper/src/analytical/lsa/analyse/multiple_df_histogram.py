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
n_param_sets=2000


turing_robustness = {}
# for nsr in [0.01, 0.05, 0.1, 0.2 ]:
for nsr in [0.01]:
    df= pickle.load( open(modellingpath + f'/3954/paper/out/analytical/lsa_dataframes/all_dataframes/lsa_df_circuit14_variantfitted7_gaussian4187715_nsr{nsr}_{n_param_sets}parametersets.pkl', "rb"))
    turing_df = df.loc[df['system_class']=='turing I oscillatory']
    turing_robustness[nsr]=(len(turing_df)/len(df))*100

# %%
import seaborn as sns
plt.figure(figsize=(4,4))
sns.barplot(x=list(turing_robustness.keys()), y=list(turing_robustness.values()), color='darkseagreen')
plt.xlabel('Relative uncertainty')
plt.ylabel('Turing I probability')
plt.title('Turing I robustness with respect to uncertainty \n Solutions around Turing fit')
plt.tight_layout()
# plt.savefig(modellingpath + f'/3954/paper/out/analytical/pyPlots/robustness_vs_uncertainty.pdf')
plt.show()

df= pickle.load( open(modellingpath + f'/3954/paper/out/analytical/lsa_dataframes/all_dataframes/lsa_df_circuit14_variantfitted7_gaussian4187715_nsr{0.01}_{n_param_sets}parametersets.pkl', "rb"))
counts = df['system_class'].value_counts()
plt.figure(figsize=(4,4))
palette = sns.color_palette("husl", n_colors=len(counts))
sns.barplot(x=counts.index, y=counts.values / len(df)*100,palette=palette)
plt.xticks(rotation=45)
plt.xlabel('System Class')
plt.ylabel('Probability')
plt.title('Frequency of solutions with 0.01 relative uncertainty around Turing I')
plt.tight_layout()
# plt.savefig(modellingpath + f'/3954/paper/out/analytical/pyPlots/robustness_vs_uncertainty.pdf')
plt.show()




#%%
# from tqdm import tqdm
# def weightTuring_single_parID(dfLoc):
#     turing_lenght = len(dfLoc.loc[dfLoc['system_class'].isin(['turing I oscillatory','turing I', 'turing I hopf'])])
#     return turing_lenght/ len(dfLoc)

# def weightTuring_full_df(df):
#     for i in tqdm(range(n_param_sets)):
#         # print(df.loc[i]['system_class'])
#         dfLoc = df.loc[i]
#         weight = weightTuring_single_parID(dfLoc)
#         df.at[i,'weightTuring']=weight
#     return df



# turing_robustness = {}
# for nsr in [0.01, 0.05, 0.1, 0.2 ]:
#     df= pickle.load( open(modellingpath + f'/3954/paper/out/analytical/lsa_dataframes/all_dataframes/lsa_df_circuit14_variantfitted7_gaussian4187715_nsr{nsr}_{n_param_sets}parametersets.pkl', "rb"))
#     # print(df['system_class'].value_counts())
#     turing_df = df.loc[df['system_class']=='turing I oscillatory']
#     turing_robustness[nsr]=(len(turing_df)/len(df))
#     weighted_df = weightTuring_full_df(df)
#     print('Turing robustness',np.sum(df['weightTuring'])/n_param_sets)


# %%
