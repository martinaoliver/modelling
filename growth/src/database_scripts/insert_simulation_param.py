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




L=50; dx =0.1; J = int(L/dx)
T =2000; dt = 0.02; N = int(T/dt)
T =10; dt = 0.1; N = int(T/dt)




#solver parameters
L=50; dx =0.05; J = int(L/dx)
T =2000; dt = 0.005; N = int(T/dt)
T =10; dt = 0.1; N = int(T/dt)



boundaryCoeff=2;rate=L/T
suggesteddt = float(dx*dx*2)
mechanism = 'openboundary'
simulation_param_dict = {'L':L, 'dx':dx, 'J':J, 'T':T, 'dt':dt, 'N':N, 
            'boundaryCoeff':boundaryCoeff, 
            'mechanism':mechanism, 'growth_rate': rate}



        
insert_simulationParam_to_sql(simulation_param_dict)