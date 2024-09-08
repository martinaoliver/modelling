

#%%
import sys
import os
pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
sys.path.append(modellingpath + '/lib')
from database.databaseFunctions import *
import pickle
#############
#############
###paths#####
#############
from tqdm import tqdm
import sys
import os
import pickle
import psycopg2
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

from numerical.cn_plot import plot1D, surfpattern


#%%
import psycopg2
credentials=f"postgresql://moliver:moliver@ld-rendres07.bc.ic.ac.uk/moliver"
with psycopg2.connect(credentials) as conn:
    with conn.cursor() as cursor:
        query = '''select ao.ss_n, mp."parID", ao."ssID", mp.variant,simulation_param_uuid,  ao.model_param_id,ao.system_class,pco.pattern_class_nogrowth from pattern_class_output pco
join analytical_output ao on (pco.model_param_id, pco."ssID") = (ao.model_param_id, ao."ssID")
join model_param mp on mp.model_param_id = ao.model_param_id
-- where (simulation_param_uuid = 'e9956efa-c4f1-439c-8146-5d009bd107d8'
where (simulation_param_uuid = '132323a4-3f93-4287-aca9-d18e84848e37'

and ( mp.variant='11' or mp.variant='12')
and mp.n_samples=1000000
and ss_n=1)

-- or (simulation_param_uuid = 'e9956efa-c4f1-439c-8146-5d009bd107d8'
or (simulation_param_uuid = '132323a4-3f93-4287-aca9-d18e84848e37'

and mp.variant='0' 
and mp.n_samples=2000000
and ss_n=1);
'''


        lsa_vs_numerical_df = pd.read_sql_query(query, conn)

lsa_vs_numerical_df

#%%

L=100; dx =0.2; J = int(L/dx)
T =18000; dt = 0.05; N = int(T/dt)
T =5000; dt = 0.05; N = int(T/dt)


L=25; dx =0.05; J = int(L/dx)
T =2000; dt = 0.005; N = int(T/dt)
T =1000; dt = 0.005; N = int(T/dt)


nested_arrays = []
n_samples=9087
for model_param_id, simulation_param_uuid in tqdm(zip(lsa_vs_numerical_df['model_param_id'].iloc[:n_samples], lsa_vs_numerical_df['simulation_param_uuid'].iloc[:n_samples])):
    print(model_param_id, simulation_param_uuid)
    U_record_1D = query_simulationOutput_single_from_sql_from_id(model_param_id,simulation_param_uuid,'U_record_1D', ssID=0)
    # cropped_U_record_1D = U_record_1D[:,-500:].astype('float16')
    cropped_U_record_1D = U_record_1D[:,-100:].astype('float16')
    if np.isnan(cropped_U_record_1D).any():
        print('nan')
        continue
    else:
        # surfpattern(cropped_U_record_1D,L,dx,J,T, morphogen=0,savefig=False,cmap='coolwarm')
        nested_arrays.append(cropped_U_record_1D)

nested_array = np.array(nested_arrays)
np.save(f'dataset_oldsim_2x500x500_crop{n_samples}.npy', nested_array)


# %%
