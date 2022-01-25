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
#execution parameters
circuit_n=2
variant=0
parametersets_n = 1000000
save_figure = False
tqdm_disable = False #disable tqdm
n_species=6
# open parameter dictionaries
general_df = pickle.load(open(modelling_home + '/3954/parameter_space_search/results/output_dataframes/lsa_df_circuit%r_variant%r_%rparametersets.pkl'%(circuit_n,variant,parametersets_n), "rb"))
# general_df = pickle.load(open(modelling_home + '/3954/parameter_space_search/results/output_dataframes/interesting_variations/ATC_5716.pkl', "rb"))
series = general_df.loc[5716]
# series = general_df.iloc[0]
# # # print(series)
df = pd.DataFrame()
kce_list = [0.01,  0.1,0.3,0.5,0.7,1,3,5,7,10,   100,1000,10000]
for i in range(len(kce_list)):
    kce = kce_list[i]
    series['kce']=kce
    df = df.append(series)
# df = df.set_index(pd.Index([1, 2, 3, 4,5,6,7]))
df = df.set_index(pd.Index(np.linspace(1,len(kce_list),len(kce_list)).astype(int)))
print(df)
print(df['kce'])
df.to_pickle('interesting_variations/ATC_5716.pkl')

# # df.loc[1,'kce']= 20
# print(df['kce'])
