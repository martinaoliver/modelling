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
import psycopg2
#%%






L=50; dx =0.1; J = int(L/dx)
T =2000; dt = 0.02; N = int(T/dt)
boundaryCoeff=1;rate=L/T
suggesteddt = float(dx*dx*2)
growth = 'nogrowth'
simulation_param_dict = {'L':L, 'dx':dx, 'J':J, 'T':T, 'dt':dt, 'N':N, 
            'boundaryCoeff':boundaryCoeff, 
            'growth':growth, 'growth rate': rate}

parID = 1544038
circuit_n='turinghill'
variant= '9'
n_samples=2000000
ssID = 0
folder = f'{circuit_n}_variant{variant}'
model_param_dict = {'parID':parID, 'circuit_n':circuit_n,'variant':variant, 'n_samples':n_samples}
filename= lambda parID: 'circuit%s_variant%s_bc%s_%s_rate%s_ID%s.%s_L%r_J%r_T%r_N%r'%(circuit_n,variant,boundaryCoeff, growth,rate,parID,ssID,L,J,T,N)
U_final_1D = pickle.load( open(modellingpath + '/growth/out/numerical/%s/simulation/%s/2Dfinal_%s.pkl'%(growth,folder,filename(parID)), 'rb'))
U_record_1D = pickle.load( open(modellingpath + '/growth/out/numerical/%s/simulation/%s/2Drecord_%s.pkl'%(growth,folder,filename(parID)), 'rb'))

#%%
print(simulation_param_dict)
query = simulationOutput_to_sql(simulation_param_dict, model_param_dict,U_final_1D,U_record_1D, ssID=ssID)
# %%
