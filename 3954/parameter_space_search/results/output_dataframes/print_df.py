import numpy
import pandas as pd
import pickle
circuit=2
variant=0
n_parametersets = 1000
df = pd.read_pickle("lsa_df_circuit%r_variant%r_%rparametersets.pkl"%(circuit,variant,n_parametersets))
pickle.dump( df.loc[0:40000], open( "lsa_df_circuit%r_variant%r_%rparametersets.pkl"%(circuit,variant,40000), "wb" ) )
