import numpy
import pandas as pd
import pickle
circuit=2
variant=0
n_parametersets = 10000
df = pd.read_pickle("lsa_df_circuit%r_variant%r_%rparametersets.pkl"%(circuit,variant,n_parametersets))
# df = pd.read_pickle('lsa_df_circuit2_variant1_1000448parametersets_rbslibrary0_concat.pkl')
# pickle.dump( df.loc[0:20000], open( "lsa_df_circuit%r_variant%r_%rparametersets.pkl"%(circuit,variant,20000), "wb" ) )
print(df.head())




print('g')
df = df.loc[8521]
print(df)        


1e-04