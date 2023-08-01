#%%
import sys
import os
pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
sys.path.append(modellingpath + '/lib')
modellingephemeral = '/rds/general/user/mo2016/ephemeral/Documents/modelling'

from database.databaseFunctions import *
import pickle
#############
#############
###paths#####
#############

import sys
import os
import pickle
import psycopg2
from tqdm import tqdm
#%%






L=50; dx =0.1; J = int(L/dx)
T =2000; dt = 0.02; N = int(T/dt)
boundaryCoeff=1;rate=L/T
suggesteddt = float(dx*dx*2)
mechanism = 'nogrowth'
simulation_param_dict = {'L':L, 'dx':dx, 'J':J, 'T':T, 'dt':dt, 'N':N, 
            'boundaryCoeff':boundaryCoeff, 
            'mechanism':mechanism, 'growth_rate': rate}

parID = 14414
circuit_n='turinghill'
variant= '9'
n_samples=2000000
ssID = 0
folder = f'{circuit_n}_variant{variant}'
model_param_dict = {'parID':parID, 'circuit_n':circuit_n,'variant':variant, 'n_samples':n_samples}
filename= lambda parID: 'circuit%s_variant%s_bc%s_%s_rate%s_ID%s_L%r_J%r_T%r_N%r'%(circuit_n,variant,boundaryCoeff, mechanism,rate,parID,L,J,T,N)
print(filename(parID), 'filename')
# %%


data_path = modellingephemeral + f'/growth/out/numerical/{mechanism}/simulation/{folder}'

parID_list = pickle.load( open(data_path + '/parID_list_%s.pkl'%(filename('x')), "rb" ) )
#%%
for parID in tqdm(parID_list):
    model_param_dict = {'parID':parID.split('.')[0], 'circuit_n':circuit_n,'variant':variant, 'n_samples':n_samples}
    U_final_1D = pickle.load( open(modellingephemeral + f'/growth/out/numerical/{mechanism}/simulation/{folder}/2Dfinal_{filename(parID)}.pkl', 'rb'))
    U_record_1D = pickle.load( open(modellingephemeral + f'/growth/out/numerical/{mechanism}/simulation/{folder}/2Drecord_{filename(parID)}.pkl', 'rb'))
    U_final_1D_list = np.array(U_final_1D).tolist()
    U_record_1D_list = np.array(U_record_1D).tolist()
    print(np.sum(U_final_1D_list))


    query = insert_simulationOutput_to_sql(simulation_param_dict, model_param_dict,U_final_1D_list,U_record_1D_list, ssID=ssID, dimensions='1D')
    # %%

