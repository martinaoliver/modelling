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

from numerical.cn_plot import plot1D, surfpattern


#%%




L=50; dx =0.1; J = int(L/dx)
T =2000; dt = 0.02; N = int(T/dt)
boundaryCoeff=1;rate=L/T
suggesteddt = float(dx*dx*2)
mechanism = 'nogrowth'
simulation_param_dict = {'L':L, 'dx':dx, 'J':J, 'T':T, 'dt':dt, 'N':N, 
            'boundaryCoeff':boundaryCoeff, 
            'mechanism':mechanism, 'growth_rate': rate}
print(simulation_param_dict.keys(), simulation_param_dict.values())
parID = 104782
circuit_n='turinghill'
variant= '9'
n_samples=2000000
ssID = 0
folder = f'{circuit_n}_variant{variant}'
model_param_dict = {'parID':parID, 'circuit_n':circuit_n,'variant':variant, 'n_samples':n_samples}



U_final_1D = query_simulationOutput_multiple_from_sql(simulation_param_dict,model_param_dict,'U_final_1D', ssID=0,fetch=2)
U_final_1D = query_simulationOutput_multiple_from_sql(simulation_param_dict,model_param_dict,'U_record_1D', ssID=0,fetch=2)
