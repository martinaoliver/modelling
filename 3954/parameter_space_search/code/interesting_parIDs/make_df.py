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
# # open parameter dictionaries
# general_df = pickle.load(open(modelling_home + '/3954/parameter_space_search/results/output_dataframes/lsa_df_circuit%r_variant%r_%rparametersets.pkl'%(circuit_n,variant,parametersets_n), "rb"))
general_df = pickle.load(open('interesting_parIDs.pkl', "rb"))
# general_df = pickle.load(open('interesting_parIDs.pkl', "rb"))
# print(general_df)
# interesting_parIDs=[16267, 14833, 20319, 25780, 29513 ,783, 20580, 2852, 8999, 25418, 1534, 15750, 19224, 25681, 6824 ,11992 ,17803, 24240, 13702, 27915, 20807, 15692, 4781, 9740 ,28593 ,5716 ,15954 ,20574 ,28641]
bullseye = [5716,25418, 17803, 28641, 2852, 28593]
# # interesting_parIDs = [1,3,4]
# general_df = pickle.load(open(modelling_home + '/3954/parameter_space_search/results/output_dataframes/interesting_variations/ATC_2852.pkl', "rb"))
# # series = general_df.loc[2852]
# bullseye_df = general_df.loc[25418]
df = general_df.loc[bullseye]
print(df)
# df = df.xs(0, level=1)# # # print(series)
#
# print(df)
# print(df['kce'])
# bullseye_df.to_pickle('bullseye_df.pkl')
df.to_pickle('bullseye_df.pkl')
print(df)
# print(bullseye_df)
# print(general_df)
# print(general_df['system_class'].value_counts())

# # df.loc[1,'kce']= 20
# print(df['kce'])
