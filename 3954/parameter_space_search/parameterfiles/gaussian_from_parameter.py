import sys
import os
pwd = os.getcwd()
root = pwd.split("home", 1)[0]
modelling_home = root + 'home/Documents/modelling'
modelling_ephemeral = root + 'ephemeral/Documents/modelling'
modulepath = modelling_home + '/3954/modules'
sys.path.append(modulepath)

from numerical_solvers_variableboundary import *
import pickle
import pandas as pd
from numpy import random
#execution parameters
circuit_n=2
variant=0
parametersets_n = 1000000
save_figure = False
tqdm_disable = False #disable tqdm
n_species=6
# open parameter dictionaries
# general_df = pickle.load(open(modelling_home + '/3954/parameter_space_search/results/output_dataframes/lsa_df_circuit%r_variant%r_%rparametersets.pkl'%(circuit_n,variant,parametersets_n), "rb"))
general_df = pickle.load(open(modelling_home + '/3954/parameter_space_search/results/output_dataframes/interesting_variations/ATC_5716.pkl', "rb"))
series = general_df.loc[11]
df = pd.DataFrame()
series = df.append(series)
df = pd.concat([series]*30000, ignore_index=True)
df = pd.concat([series]*100, ignore_index=True)

#ADD NORMAL DISTRIBUTIONS
parameter_list = ['Va','Vb','Vc','Vd','Ve','Vf','ba','bb','bc','bd','be','bf','d_B','kaa','kbd','kce','kda','keb','kee','kfe','mua','mulva']
for parameter in parameter_list:
    df[parameter] = random.normal(loc=df[parameter].iloc[0], scale=df[parameter].iloc[0]*0.25, size=len(df))

 # [234] * len(df)
del df['ss_class']
del df['ss_n']
del df['system_class']
del df['ss_list']
del df['maxeig']
# print(df)
print(df)
# df.to_pickle('df_circuit2_variant5716gaussian_30000parametersets.pkl')
