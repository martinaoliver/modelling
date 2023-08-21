#%%


import sys
import os
pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
modellingephemeral = '/rds/general/ephemeral/user/mo2016/ephemeral/Documents/modelling'
sys.path.append(modellingpath + '/lib')

from database.databaseFunctions import *

import pickle
#############
#############
###paths#####
#############

import sys
import os
import pickle
import numpy as np
from tqdm import tqdm
#%%





# slow
L=20; dx =0.1; J = int(L/dx)
T =100; dt = 0.02; N = int(T/dt)
boundaryCoeff = 1
division_time_hours=0.5
p_division=0.38;seed=1


# # medium
# L=20; dx =0.1; J = int(L/dx)
# T =50; dt = 0.02; N = int(T/dt)
# boundaryCoeff = 1
# division_time_hours=0.5
# p_division=1;seed=1



# # fast
L=20; dx =0.1; J = int(L/dx)
T =25; dt = 0.02; N = int(T/dt)
boundaryCoeff = 1
division_time_hours=0.2
p_division=0.7;seed=1

#caFastMotherMultiple
L=25; dx =0.1; J = int(L/dx)
T =110; dt = 0.02; N = int(T/dt)
boundaryCoeff = 1
division_time_hours=0.5
p_division=0.38;seed=1

shape='caFastMotherMultiple'
simulation_param_dict = {'L':L, 'dx':dx, 'J':J, 'T':T, 'dt':dt, 'N':N, 
            'boundaryCoeff':boundaryCoeff, 
            'shape':shape, 'p_division': p_division, 'seed':seed, 'division_time_hours':division_time_hours}




ssID=0
circuit_n=14 #circuit_n='circuit14'
variant='fitted7_gaussian4187715_nsr0.01' #'2nd'#'fitted7_gaussian4187715_nsr0.01' #variant='fitted7'
Kce=100
n_samples = 2000 #1000000#2000 #n_samples = 13700000
folder = 'circuit14variantfitted7_gaussian4187715'#'circuit14variant2ndBalancedKce100'#'circuit14variantfitted7_gaussian4187715'


#ssID=0
#circuit_n=14 #circuit_n='circuit14'
#variant='fitted7_gaussian4187715_nsr0.01' #variant='fitted7'
#Kce=100
#n_samples =2000 #n_samples = 13700000
#folder ='circuit14variantfitted7_gaussian4187715'


model_param_dict = {'circuit_n':circuit_n,'variant':variant, 'n_samples':n_samples}
filename= lambda parID: 'circuit%r_variant%s_bc%s_%s_ID%s_L%r_J%r_T%r_N%r'%(circuit_n,variant,boundaryCoeff, shape,parID,L,J,T,N)
# filename = lambda parID: f'circuit{circuit_n}_variant{variant}_{n_samples}parametersets_balanced_Kce{Kce}_bc{boundaryCoeff}_{shape}_ID{parID}_L{L}_J{J}_T{T}_N{N}'



# ssID=0
# circuit_n='14' #circuit_n='circuit14'
# variant='2nd' #variant='fitted7'
# balance='balanced'
# Kce=100
# n_samples = 1000000 #n_samples = 13700000
# folder = 'circuit14variant2ndBalancedKce100'
# filename= lambda parID: f'circuit{circuit_n}_variant{variant}_{n_samples}parametersets_{balance}_Kce{Kce}_bc{boundaryCoeff}_{shape}_ID{parID}_L{L}_J{J}_T{T}_N{N}'
# filename= lambda parID: f'circuit{circuit_n}_variant{variant}_{n_samples}parametersets_{balance}_Kce{Kce}_bc{boundaryCoeff}_{shape}_ID{parID}_L{L}_J{J}_T{T}_N{N}'



data_path = modellingpath + '/3954/paper/out/numerical/colonies/simulation/%s'%(folder)
parID_list = pickle.load( open(data_path + '/parID_list_%s.pkl'%(filename('x')), "rb" ) )
#%%
parID_list = [141318]
for parID in tqdm(parID_list[:2]):
    print(parID)
    
    model_param_dict = {'parID':parID, 'circuit_n':circuit_n,'variant':variant, 'n_samples':n_samples}
    print('filename', filename(parID))
    U_final_1D = pickle.load( open(modellingpath + f'/3954/paper/out/numerical/colonies/simulation/{folder}/2Dfinal_{filename(parID)}.pkl', 'rb'))
    U_record_1D = pickle.load( open(modellingpath + f'/3954/paper/out/numerical/colonies/simulation/{folder}/2Drecord_{filename(parID)}.pkl', 'rb'))
    U_final_1D_list = np.array(U_final_1D).tolist()
    U_record_1D_list = np.array(U_record_1D).tolist()
    print(np.sum(U_final_1D_list))


    query = insert_simulationOutput_to_sql(simulation_param_dict, model_param_dict,U_final_1D_list,U_record_1D_list, ssID=ssID, dimensions='1D')



# %%
parID=6
model_param_dict = {'parID':parID, 'circuit_n':circuit_n,'variant':variant, 'n_samples':n_samples}
print('filename', filename(parID))
U_final_1D = pickle.load( open(modellingpath + f'/3954/paper/out/numerical/colonies/simulation/{folder}/2Dfinal_{filename(parID)}.pkl', 'rb'))
U_record_1D = pickle.load( open(modellingpath + f'/3954/paper/out/numerical/colonies/simulation/{folder}/2Drecord_{filename(parID)}.pkl', 'rb'))
U_final_1D_list = np.array(U_final_1D).tolist()
U_record_1D_list = np.array(U_record_1D).tolist()
print(np.sum(U_final_1D_list))


query = insert_simulationOutput_to_sql(simulation_param_dict, model_param_dict,U_final_1D_list,U_record_1D_list, ssID=ssID, dimensions='1D')


# %%
