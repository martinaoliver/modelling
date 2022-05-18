import numpy
import pandas as pd
import pickle
circuit=2
variant=0
n_parametersets = 1000000
n_parametersets = 10000
# df = pd.read_pickle("df_circuit%r_variant%r_%rparametersets.pkl"%(circuit,variant,n_parametersets))
df = pd.read_pickle('df_circuit2_variant1_1954parametersets_rbslibrary0.pkl')

print('g')
# pickle.dump( df.loc[0:20000], open( "df_circuit%r_variant%r_%rparametersets.pkl"%(circuit,variant,10000), "wb" ) )
# print(df.loc[5716])
print(df)