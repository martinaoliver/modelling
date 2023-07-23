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
import matplotlib.pyplot as plt

from numerical.plotting_numerical import *


#%%

# # slow
# L=20; dx =0.1; J = int(L/dx)
# T =100; dt = 0.02; N = int(T/dt)
# boundaryCoeff = 1
# divisionTimeHours=0.5
# p_division=0.38;seed=1


# # medium
# L=20; dx =0.1; J = int(L/dx)
# T =50; dt = 0.02; N = int(T/dt)
# boundaryCoeff = 1
# divisionTimeHours=0.5
# p_division=1;seed=1



# fast
L=20; dx =0.1; J = int(L/dx)
T =25; dt = 0.02; N = int(T/dt)
boundaryCoeff = 1
divisionTimeHours=0.2
p_division=0.7;seed=1
x_gridpoints=int(1/dx)


shape='ca'
simulation_param_dict = {'L':L, 'dx':dx, 'J':J, 'T':T, 'dt':dt, 'N':N, 
            'boundaryCoeff':boundaryCoeff, 
            'shape':'ca', 'p_division': p_division, 'seed':seed}





ssID=0
circuit_n='14' #circuit_n='circuit14'
variant='2nd' #variant='fitted7'
balance='Balanced'
Kce=100
n_samples = 1000000 #n_samples = 13700000
folder = 'circuit14variant2ndBalancedKce100'
model_param_dict = {'parID':1, 'circuit_n':circuit_n,'variant':variant, 'n_samples':n_samples, 'balance':balance}



#%%



#%%
simulationOutput = query_simulationOutput_from_sql(simulation_param_dict, model_param_dict,query_column = 'U_final_1D', ssID=ssID)
plot_redgreen_contrast(simulationOutput,L)
# %%ssID=0
circuit_n='14' #circuit_n='circuit14'
variant='2nd' #variant='fitted7'
balance='Balanced'
Kce=100
n_samples = 1000000 #n_samples = 13700000
folder = 'circuit14variant2ndBalancedKce100'
model_param_dict = {'parID':1, 'circuit_n':circuit_n,'variant':variant, 'n_samples':n_samples, 'balance':balance}



#%%



#%%
simulationOutput = query_simulationOutput_from_sql(simulation_param_dict, model_param_dict,query_column = 'U_final_1D', ssID=ssID)
plot_redgreen_contrast(simulationOutput,L)
# %%