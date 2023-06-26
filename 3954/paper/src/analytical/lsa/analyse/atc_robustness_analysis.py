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
import seaborn as sns
from scipy.optimize import curve_fit

#%%

# Specify name of circuit and variant investigated
circuit_n='circuit14'
variant='2nd'
balance = 'balanced'
# Specifiy number of parameter sets in parameterset file to be loaded
n_param_sets = 1000000

print(f'{circuit_n}, Variant:{variant}, {balance}')
Kce_list = [0.1,1,10,100,1000]
turingStates = ['turing I','turing I oscillatory', 'turing I hopf']  
Kce_robustness_dict = {}
for Kce in Kce_list:
    df= pickle.load( open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/all_dataframes/lsa_df_%s_variant%s_%rparametersets_%s_Kce%s.pkl'%(circuit_n,variant,n_param_sets,balance,Kce), "rb"))
    turing_df = df.loc[df['system_class'].isin(turingStates)]
    robustness_percentage = len(turing_df)/len(df)*100
    Kce_robustness_dict[Kce] = robustness_percentage
    print(Kce, robustness_percentage)
# %%
robustness_percentage_list = list(Kce_robustness_dict.values())

# plt.scatter(Kce_list, robustness_percentage_list)
def objective(x, a, b, c):
    return a * np.exp(-b * x) + c

popt, _ = curve_fit(objective, Kce_list, robustness_percentage_list)
a, b, c = popt
# use optimal parameters to calculate new values
x_new= np.logspace(-1,3,100)
y_new = objective(x_new, a, b, c)
plt.plot(x_new, y_new)
plt.fill_between(x_new, y_new,  alpha=0.4)

plt.xscale('log')
plt.xlabel('Kce parameter')
plt.ylabel('% of Turing robustness')
plt.tight_layout()
plt.savefig(modellingpath + '/3954/paper/out/analytical/pyPlots/Kce_robustness_df_%s_variant%s_%rparametersets.pdf'%(circuit_n,variant,n_param_sets))
plt.savefig(modellingpath + '/3954/paper/out/analytical/pyPlots/Kce_robustness_df_%s_variant%s_%rparametersets.png'%(circuit_n,variant,n_param_sets))
plt.show()
# %%









#%%

# Specify name of circuit and variant investigated
circuit_n='circuit14'
variant='2nd'
balance = 'balanced'
# Specifiy number of parameter sets in parameterset file to be loaded
n_param_sets = 1000000

print(f'{circuit_n}, Variant:{variant}, {balance}')

# df= pickle.load( open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/all_dataframes/lsa_df_%s_variant%s_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), "rb"))
# df= pickle.load( open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/all_dataframes/lsa_df_%s_variant%s_%rparametersets_balancedSemiBalanced.pkl'%(circuit_n,variant,n_param_sets), "rb"))
df= pickle.load( open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/all_dataframes/lsa_df_%s_variant%s_%rparametersets_%s.pkl'%(circuit_n,variant,n_param_sets,balance), "rb"))
# df= pickle.load( open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/all_dataframes/lsa_df_%s_variant%s_%rparametersets_balancedSemiBalanced.pkl'%(circuit_n,variant,n_param_sets), 'rb'))

print(df['system_class'].value_counts())
#%%

sns.kdeplot(list(df.loc[:,'Kce'].values), fill=True,log_scale=True,cut=1,color='teal', linewidth = 4, alpha = 0.1,label = 'HSL Diffusion')

