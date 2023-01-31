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
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from equations.class_subcircuit_eq import *
from sklearn.metrics import mean_squared_error
from tqdm import tqdm


def checkBalance(par_dict):
    balanceDict = {}
    for Km in Km_list:
        # print(Km)
        Vx =par_dict[KtoV[Km]]
        Kxy = par_dict[Km]
        if Kxy >= 1 and Kxy <= Vx:
            balanceDict[Km] = 'Balanced'
        elif Kxy > 0.1 and Kxy < Vx*10:
            balanceDict[Km] ='Semi balanced'
        elif Kxy <= 0.1 or Kxy >= Vx*10:
            balanceDict[Km] ='Not balanced'
        else:
            print('ERROR!!!!!!!!!')

    if 'Not balanced' in balanceDict.values():
        return 'Not balanced'
    elif 'Semi balanced'  in balanceDict.values():
        return 'Semi balanced'
    elif all(x == 'Balanced' for x in balanceDict.values()):
        return 'Balanced'




# Specify name of circuit and variant investigated
circuit_n='circuit14'
variant='1nd'
# Specifiy number of parameter sets in parameterset file to be loaded
n_param_sets = 5000000

print(f'Circuit:{circuit_n}, Variant:{variant}')

df_full= pkl.load( open(modellingpath + "/3954/paper/input/lhs_parameterfiles/df_%s_variant%s_%rparametersets.pkl"%(circuit_n,variant,n_param_sets), "rb"))





Km_list = ['Kda', 'Kab', 'Keb', 'Kbd', 'Kfe', 'Kee', 'Kce' ]
KtoV = {'Kda': 'VD', 'Kab': 'VA', 'Keb': 'VE', 'Kbd': 'VB', 'Kfe': 'VF', 'Kee': 'VE', 'Kce': 'VC' }
balanceList = []    
for parID in tqdm(df_full.index):
    par_dict = df_full.loc[parID].to_dict()
    balanceList.append(checkBalance(par_dict))
df_full['balance'] = balanceList



print(df_full['balance'].value_counts())
# df_full[df_full['balance'] == 'Balanced']
# pkl.dump(df_full, open(modellingpath + "/3954/paper/input/balanced_parameterfiles/df_%s_variant%s_%rparametersets_balanced.pkl"%(circuit_n,variant,n_param_sets), "wb"))
