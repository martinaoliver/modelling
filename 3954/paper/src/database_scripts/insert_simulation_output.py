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

import sys
import os
import pickle
import numpy as np
from tqdm import tqdm
#%%





# # slow
# L=20; dx =0.1; J = int(L/dx)
# T =100; dt = 0.02; N = int(T/dt)
# boundaryCoeff = 1
# divisionTimeHours=0.5
# p_division=0.38;seed=1


# medium
L=20; dx =0.1; J = int(L/dx)
T =50; dt = 0.02; N = int(T/dt)
boundaryCoeff = 1
divisionTimeHours=0.5
p_division=1;seed=1



# # fast
# L=20; dx =0.1; J = int(L/dx)
# T =25; dt = 0.02; N = int(T/dt)
# boundaryCoeff = 1
# divisionTimeHours=0.2
# p_division=0.7;seed=1


shape='ca'
simulation_param_dict = {'L':L, 'dx':dx, 'J':J, 'T':T, 'dt':dt, 'N':N, 
            'boundaryCoeff':boundaryCoeff, 
            'shape':'ca', 'p_division': p_division, 'seed':seed}



ssID=0
circuit_n='14' #circuit_n='circuit14'
variant='2nd' #variant='fitted7'
balance='balanced'
Kce=100
n_samples = 1000000 #n_samples = 13700000
folder = 'circuit14variant2ndBalancedKce100'
filename= lambda parID: f'circuit{circuit_n}_variant{variant}_{n_samples}parametersets_{balance}_Kce{Kce}_bc{boundaryCoeff}_{shape}_ID{parID}_L{L}_J{J}_T{T}_N{N}'

data_path = modellingpath + '/3954/paper/out/numerical/colonies/simulation/%s'%(folder)
parID_list = pickle.load( open(data_path + '/parID_list_%s.pkl'%(filename('x')), "rb" ) )
#%%
for parID in parID_list:
    model_param_dict = {'parID':parID, 'circuit_n':circuit_n,'variant':variant, 'n_samples':n_samples, 'balance':balance}
    U_final_1D = pickle.load( open(modellingpath + f'/3954/paper/out/numerical/colonies/simulation/{folder}/2Dfinal_{filename(parID)}.pkl', 'rb'))
    U_record_1D = pickle.load( open(modellingpath + f'/3954/paper/out/numerical/colonies/simulation/{folder}/2Drecord_{filename(parID)}.pkl', 'rb'))
    U_final_1D_list = np.array(U_final_1D).tolist()
    U_record_1D_list = np.array(U_record_1D).tolist()
    #%%

    query = simulationOutput_to_sql(simulation_param_dict, model_param_dict,U_final_1D_list,U_record_1D_list, ssID=ssID)
    # %%

