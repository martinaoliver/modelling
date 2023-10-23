#%%

import sys
import os

pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
sys.path.append(modellingpath + '/lib')

from database.databaseFunctions import *
import pickle
import psycopg2
import time
import numpy as np
from tqdm import tqdm

#%%
#############

# #%%
# Specify name of circuit and variant investigated



circuit_n='turinghill'
variant=int(sys.argv[1])
seed=int(sys.argv[2])
n_species=2
# Specifiy number of parameter sets in parameterset file to be loaded
n_samples = int(sys.argv[3])
threshold=int(sys.argv[4])

# circuit_n='turinghill'
# variant=11
# seed=0
# n_species=2
# # Specifiy number of parameter sets in parameterset file to be loaded
# n_samples =10
# threshold=3000



lhs_df = pickle.load( open(modellingpath + '/growth/input/parameterfiles/df_%s_variant%s_%rparametersets_seed%s.pkl'%(circuit_n,variant,n_samples,seed), "rb"))
lsa_df = pickle.load( open(modellingpath + '/growth/out/analytical/lsa_dataframes/lsa_df_%s_variant%s_%rparametersets_seed%s.pkl'%(circuit_n,variant,n_samples, seed), "rb"))

#%%
query = lambda system_class: f'''select count(*) from analytical_output
join model_param mp on mp.model_param_id = analytical_output.model_param_id
where circuit_n='{circuit_n}'
and variant='{variant}'
and n_samples={n_samples}
and system_class='{system_class}'
'''

for system_class in list(lsa_df['system_class'].unique()):
    n_system_class_ssn1 = general_query(query(system_class))
    print(system_class, n_system_class_ssn1[0][0][0])
    if  n_system_class_ssn1[0][0][0]<threshold:
        selected_lsa_df = lsa_df.loc[lsa_df['system_class']==system_class]
        selected_lsa_df = selected_lsa_df.iloc[:threshold]
        print(len(selected_lsa_df))

        selected_lhs_df= lhs_df.loc[selected_lsa_df.index.get_level_values(0)]
        selected_lhs_df = selected_lhs_df.drop_duplicates()
        
        modelParam_df_to_sql(selected_lhs_df, circuit_n, variant, n_samples)
        analyticalOutput_df_to_sql(selected_lsa_df, circuit_n, variant, n_samples)



# %%
