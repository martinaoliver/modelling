
#%%
###paths#####
#############
import sys
import os

from importlib_metadata import distribution
pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
sys.path.append(modellingpath + '/lib')
#############

from analytical.linear_stability_analysis import big_turing_analysis_df, detailed_turing_analysis_dict
from randomfunctions import plot_all_dispersion
from database.databaseFunctions import *

import pickle
import numpy as np 
#######################
#########CODE##########
#######################

circuit_n='turinghill'; n_species=2; variant = '9'; n_samples=1000000
parameters_df= pickle.load( open(modellingpath + f'/growth/input/parameterfiles/df_turinghill_variant{variant}_2000000parametersets.pkl', 'rb'))
# df= pickle.load( open(modellingpath + f'/growth/out/analytical/instability/instability_df_circuitturinghill_variant{variant}_combinedparametersets.pkl', 'rb'))
df= pickle.load( open(modellingpath + f'/growth/out/analytical/lsa_dataframes/lsa_df_circuitturinghill_variant{variant}_combinedparametersets.pkl', "rb" ) )
df = df[~df.index.duplicated(keep='first')]

# pickle.dump(result_df, open( modellingpath + f'/growth/out/analytical/lsa_dataframes/lsa_df_circuitturinghill_variant8-9_combinedparametersets.pkl', "wb" ) )

#make df like lhs_df
df = df.reindex(columns=parameters_df.columns)

# df = df.iloc[:10]

# %%
output_df = big_turing_analysis_df(df,circuit_n,n_species,print_parID=False, tqdm_disable=False)
output_df
# %%
tupled_index =  [tuple(l) for l in output_df.index]
multi_index = pd.MultiIndex.from_tuples(tupled_index)
output_df = output_df.set_index(multi_index)

#%%
pickle.dump(output_df, open( modellingpath + f'/growth/out/analytical/lsa_dataframes/multiinstability_lsa_df_circuitturinghill_variant{variant}_combinedparametersets.pkl', "wb" ) )

#%%


analyticalOutput_df_to_sql(output_df, circuit_n, variant, n_samples)


# %%
