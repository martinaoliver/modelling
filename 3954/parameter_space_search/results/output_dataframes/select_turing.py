import numpy
import pandas as pd
import pickle
circuit=2
variant=9
n_parametersets = 1000000
df = pd.read_pickle("lsa_df_circuit%r_variant%r_%rparametersets.pkl"%(circuit,variant,n_parametersets))
# df = pd.read_pickle('lsa_df_circuit2_variant1_1000448parametersets_rbslibrary0_concat.pkl')
# print(df.head())

states = ['turing I'] 
turing = df.loc[df['system_class'].isin(states)]
# print(turing.index.levels[1])
turing = turing.xs(0, level=1)
print(turing)
# print(turing.index.get_level_values(0))
print('g')
# df = df.loc[5716]
pickle.dump( turing, open( "../turing_dataframes/turing_lsa_df_circuit%r_variant%r_%rparametersets.pkl"%(circuit,variant,n_parametersets), "wb" ) )
