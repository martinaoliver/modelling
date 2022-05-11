import numpy
import pandas as pd
import pickle
circuit=2
variant=0
n_parametersets = 1000000
# df = pd.read_pickle("lsa_df_circuit%r_variant%r_%rparametersets.pkl"%(circuit,variant,n_parametersets))
df = pd.read_pickle('lsa_df_circuit2_variant0_51200parametersets_batch9_rbslibrary0.pkl')
# pickle.dump( df.loc[0:20000], open( "lsa_df_circuit%r_variant%r_%rparametersets.pkl"%(circuit,variant,20000), "wb" ) )
# print(df.head())
# print(df['system_class'].value_counts())
print('g')
# df = df.loc[5716]
print(df)