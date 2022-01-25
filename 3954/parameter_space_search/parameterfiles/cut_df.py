import pandas as pd
import pickle
print('load df')
df= pickle.load( open('df_circuit2_variant0_1000000parametersets.pkl', "rb" ) )
print('df loaded')
df_cut =df.iloc[:10000]
df_cut.to_pickle('df_circuit2_variant0_10000parametersets.pkl')
