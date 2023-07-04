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

def simulationOutput_df_to_sql(sim_dict,):
    #query simulation_param_id using sim_dict
    #query model_param_id using circuit, variant, nsamples, 
    #query steadystate. default 0 for now
    

    return



L=50; dx =0.1; J = int(L/dx)
T =2000; dt = 0.02; N = int(T/dt)
boundaryCoeff=1;rate=L/T
suggesteddt = float(dx*dx*2)

sim_dict = {'L':L, 'dx':dx, 'J':J, 'T':T, 'dt':dt, 'N':N, 
            'boundaryCoeff':boundaryCoeff, 
            'growth':'nogrowth', 'growth rate': rate}
simulationOutput_df_to_sql(sim_dict)

