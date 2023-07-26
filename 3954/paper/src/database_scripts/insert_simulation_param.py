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





# slow
L=20; dx =0.1; J = int(L/dx)
T =100; dt = 0.02; N = int(T/dt)
boundaryCoeff = 2
division_time_hours=0.5
p_division=0.38;seed=1


# # medium
# L=20; dx =0.1; J = int(L/dx)
# T =50; dt = 0.02; N = int(T/dt)
# boundaryCoeff = 2
# division_time_hours=0.5
# p_division=1;seed=1

# T= 1;dt= 0.5; N= 2,

# # fast
# L=20; dx =0.1; J = int(L/dx)
# T =25; dt = 0.02; N = int(T/dt)
# boundaryCoeff = 1
# division_time_hours=0.2
# p_division=0.7;seed=1

shape='ca'
sim_dict = {'L':L, 'dx':dx, 'J':J, 'T':T, 'dt':dt, 'N':N, 
            'boundaryCoeff':boundaryCoeff, 
            'shape':'ca', 'p_division': p_division,  'seed':seed, 'division_time_hours': division_time_hours}

insert_simulationParam_to_sql(sim_dict)




