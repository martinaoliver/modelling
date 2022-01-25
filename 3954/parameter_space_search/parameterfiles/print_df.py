import numpy
import pandas as pd
import pickle
circuit=2
variant=0
n_parametersets = 1000000
df = pd.read_pickle("df_circuit%r_variant%r_%rparametersets.pkl"%(circuit,variant,n_parametersets))
df = df.iloc[:1000]
pickle.dump(df, open('df_circuit%r_variant%r_%rparametersets.pkl'%(circuit,variant,len(df)), 'wb'))

# print(df)
