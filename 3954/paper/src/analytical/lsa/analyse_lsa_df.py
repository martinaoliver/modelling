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
import pandas as pd

# Specify name of circuit and variant investigated
circuit_n='circuit2'
variant=0
n_species=6
# Specifiy number of parameter sets in parameterset file to be loaded
n_param_sets = 1000000

print(f'Circuit:{circuit_n}, Variant:{variant}')

df= pickle.load( open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/all_dataframes/lsa_df_%s_variant%r_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), "rb"))
print(df['system_class'].value_counts())

# #values for which complex dispersion = true
# complex_df = df[df['complex_dispersion']==True]
# complex_df.index  = complex_df.index.droplevel(-1)

# #values that have instabilities
instabilities = ['turing I', 'turing II', 'turing I hopf', 'turing I oscillatory', 'turing II hopf','hopf', 'turing semi-hopf']  
instabilities_df = df.loc[df['system_class'].isin(instabilities)]
instabilities_df.index  = instabilities_df.index.droplevel(-1)
print(instabilities_df['system_class'].value_counts())
pickle.dump( instabilities_df, open(modellingpath + '/3954/paper/out/analytical/lsa_dataframes/instabilities_dataframes/instability_df_%s_variant%r_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), "wb" ) )
len(instabilities_df)
# instabilityComplex_df = pd.concat([complex_df, instabilities_df])
# pickle.dump( instabilityComplex_df, open(modellingpath + '/growth/out/analytical/instabilityComplex/instabilityComplex_df_%s_variant%r_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), "wb" ) )

# #values that have turing
# turingStates = ['turing I','turing I oscillatory']  
# turing_df = df.loc[df['system_class'].isin(turingStates)]
# # turing_df = turing_df.xs(0, level=1)
# turing_df.index  = turing_df.index.droplevel(-1)
# pickle.dump( turing_df, open(modellingpath + '/growth/out/analytical/turing/turing_df_%s_variant%r_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), "wb" ) )
# #remove secondary index from turing_df


# # select_turing=True
# # if select_turing == True:    
# #     states = ['turing I, turing II', 'turing I hopf', 'turing I oscillatory', 'turing II hopf','hopf']  
# #     turing_df = df.loc[df['system_class'].isin(states)]
# #     turing_df = turing_df.xs(0, level=1)
# #     pickle.dump( turing_df, open(modellingpath + '/growth/out/analytical/turing_dataframes/turing_df_%s_variant%r_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), "wb" ) )

# # select_turingI=True
# # if select_turingI == True:
# #     states = ['turing I', 'turing I oscillatory']
# #     turingI_df = df.loc[df['system_class'].isin(states)]
# #     print(turingI_df)
# #     if len(turingI_df) > 0:
# #         turingI_df = turingI_df.xs(0, level=1)
# #         pickle.dump( turingI_df, open(modellingpath + '/growth/out/analytical/turing_dataframes/turingI_df_%s_variant%r_%rparametersets.pkl'%(circuit_n,variant,n_param_sets), "wb" ) )
# # crop_to = 100000
# # cropped_df = df.iloc[:crop_to]
