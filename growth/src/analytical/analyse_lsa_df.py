#############
###paths#####
#############
import sys
import os

from importlib_metadata import distribution
pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
sys.path.append(modellingpath + '/lib')
#############


import pickle

circuit_n='turinghill'
variant= 1
n_species=2
# Specifiy number of parameter sets in parameterset file to be loaded
n_param_sets = 2000000


df= pickle.load( open(modellingpath + '/growth/out/analytical/lsa_dataframes/lsa_df_%s_variant%r_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), "rb"))
print(df['system_class'].value_counts())

select_turing=True
if select_turing == True:
    states = ['turing I, turing II', 'turing I hopf', 'turing I oscillatory', 'turing II hopf','hopf']  
    turing_df = df.loc[df['system_class'].isin(states)]
    turing_df = turing_df.xs(0, level=1)
    pickle.dump( turing_df, open(modellingpath + '/growth/out/analytical/turing_dataframes/turing_df_%s_variant%r_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), "wb" ) )

select_turingI=True
if select_turingI == True:
    states = ['turing I']  
    turingI_df = df.loc[df['system_class'].isin(states)]
    print(turingI_df)
    if len(turingI_df) > 0:
        turingI_df = turingI_df.xs(0, level=1)
        pickle.dump( turingI_df, open(modellingpath + '/growth/out/analytical/turing_dataframes/turingI_df_%s_variant%r_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), "wb" ) )
crop_to = 100000
cropped_df = df.iloc[:crop_to]
pickle.dump( cropped_df, open(modellingpath + '/growth/out/analytical/lsa_dataframes/lsa_df_%s_variant%r_%rparametersets.pkl'%(circuit_n,variant,100000), "wb" ) )