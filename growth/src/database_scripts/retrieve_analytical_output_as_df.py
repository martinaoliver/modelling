#%%
import sys
import os

pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
sys.path.append(modellingpath + '/lib')

from database.databaseFunctions import *
import pickle

#%%


circuit_n='\'turinghill\'';variant='\'8\'';n_species=2;n_samples=1000000
model_param_dict = {'circuit_n': circuit_n, 'variant': variant , 'n_samples':n_samples}
result_df_8 = query_analyticalOutput_df_from_sql(model_param_dict)
result_df_8

#%%

circuit_n='\'turinghill\'';variant='\'9\'';n_species=2;n_samples=1000000
model_param_dict = {'circuit_n': circuit_n, 'variant': variant , 'n_samples':n_samples}
result_df_9 = query_analyticalOutput_df_from_sql(model_param_dict)
result_df_9
# %%
pickle.dump(result_df_8, open( modellingpath + f'/growth/out/analytical/lsa_dataframes/lsa_df_circuitturinghill_variant8_combinedparametersets.pkl', "wb" ) )
pickle.dump(result_df_9, open( modellingpath + f'/growth/out/analytical/lsa_dataframes/lsa_df_circuitturinghill_variant9_combinedparametersets.pkl', "wb" ) )

#%%


result_df = pd.concat([result_df_8,result_df_9])

pickle.dump(result_df, open( modellingpath + f'/growth/out/analytical/lsa_dataframes/lsa_df_circuitturinghill_variant8-9_combinedparametersets.pkl', "wb" ) )

print(result_df.columns)

#%%

instabilities = ['turing I', 'turing II', 'turing I hopf', 'turing I oscillatory', 'turing II hopf','hopf', 'turing semi-hopf']  
result_df_instabilities = result_df.loc[result_df['system_class'].isin(instabilities)]
pickle.dump(result_df_instabilities, open( modellingpath + f'/growth/out/analytical/instability/instability_df_circuitturinghill_variant8-9_combinedparametersets.pkl', "wb" ) )

result_df_instabilities_8 = result_df_8.loc[result_df_8['system_class'].isin(instabilities)]
pickle.dump(result_df_instabilities_8, open( modellingpath + f'/growth/out/analytical/instability/instability_df_circuitturinghill_variant8_combinedparametersets.pkl', "wb" ) )

result_df_instabilities_9 = result_df_9.loc[result_df_9['system_class'].isin(instabilities)]
pickle.dump(result_df_instabilities_9, open( modellingpath + f'/growth/out/analytical/instability/instability_df_circuitturinghill_variant9_combinedparametersets.pkl', "wb" ) )


# %%
print(result_df_instabilities['system_class'].value_counts())

# %%
instabilities_df1= pickle.load( open(modellingpath + f'/growth/out/analytical/instability/instability_df_circuitturinghill_variant8-9_combinedparametersets.pkl', 'rb'))

# %%
print(instabilities_df1['system_class'].value_counts())

# %%
