import numpy
import pandas as pd
import pickle
circuit=2
variant=0
new_variant = 5
n_parametersets = 1000000
df = pd.read_pickle("df_circuit%r_variant%r_%rparametersets.pkl"%(circuit,variant,n_parametersets))
df['ba'],df['bb'],df['bc'],df['bd'],df['be'],df['bf']=[0.5]*6
# df['n']=3
pickle.dump(df, open('df_circuit%r_variant%r_%rparametersets.pkl'%(circuit,new_variant,len(df)), 'wb'))
print(df.iloc[:10])

# df['ba']
