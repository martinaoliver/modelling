#%%#############
###paths#####
#############
import sys
import os

pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
sys.path.append(modellingpath + '/lib')
#############

from numerical.cn_edgegrowth2_numba import cn_edgegrowth2 as cn_edgegrowth2_numba
from numerical.cn_edgegrowth2 import cn_edgegrowth2
from numerical.cn_nogrowth import cn_nogrowth
from numerical.cn_nogrowth_numba import cn_nogrowth_numba

from numerical.cn_plot import plot1D, surfpattern
from database.databaseFunctions import *

import pickle
import matplotlib.pyplot as plt
import time
import numpy as np

#system parameters
circuit_n = 'turinghill'
variant=9
n_samples = 2000000

# par_dicxt = {'c1':0.1, 'c2':1,'c3':0.9,'c4':1, 'd_A': 1, 'd_B':10}
df= pickle.load( open(modellingpath + "/growth/out/analytical/turing/turing_df_%s_variant%r_%rparametersets.pkl"%(circuit_n,variant,n_samples), "rb"))

# df = multiple_df.xs(0, level=1)
#solver parameters
L=10; dx =0.1; J = int(L/dx)
T =20; dt = 0.02; N = int(T/dt)

# L=10; dx =1; J = int(L/dx)
# T =30; dt = 0.5; N = int(T/dt)

# L=50; dx =0.1; J = int(L/dx)
# T =2000; dt = 0.02; N = int(T/dt)
boundaryCoeff=1;rate=L/T
suggesteddt = float(dx*dx*2)



suggesteddt = float(dx*dx*2)

print(f'suggested dt = {suggesteddt}, used dt = {dt}')

suggesteddt = float(dx*dx*2)
print(f'suggested dt = {suggesteddt}, used dt = {dt}')

# for parID,ss in df.index:
# parID= (14414,0) #parameter set to use
parID=1930331 ;ssID=2#parameter set to use
par_dict = df.loc[parID,ssID].to_dict()
# ssID=par_dict['']

print(par_dict)
parameter_to_modify = ['ba', 'bb', 'Va', 'Vb', 'mua', 'mub']
for parameter in parameter_to_modify:
    par_dict[parameter] = par_dict[parameter] 
print(par_dict)
# par_dict = df.loc[parID].to_dict()
print(f'estimated wavelenght = {par_dict["estimated_wvl"]}')


print(par_dict)
model_param_dict = {'parID':parID, 'circuit_n':circuit_n,'variant':variant, 'n_samples':n_samples}



#%%
#run

growth = True
if growth == True:
    growth = 'edgegrowth3'
    boundaryCoeff=2
    simulation_param_dict = {'L':L, 'dx':dx, 'J':J, 'T':T, 'dt':dt, 'N':N, 'boundaryCoeff':boundaryCoeff, 'growth':growth, 'growth rate': rate}
    st = time.time()
    U_final_1D,U_record_1D, U0, x_grid, reduced_t_grid, cellMatrix= cn_edgegrowth3(par_dict,L,J,T,N, circuit_n ,rate=rate, boundaryCoeff=boundaryCoeff, tqdm_disable=True)
    elapsed_time = time.time() - st
    print('Execution time no numba:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
    plt.scatter(x_grid,cellMatrix)
    plt.show()
    

    #plot
    plot1D(U_final_1D, savefig=False,filename='')
    plt.show()
    surfpattern(U_record_1D, [x_grid, reduced_t_grid], 'linear',morphogen=0, rate=0, savefig=False,filename='',logResults=False,normalize=False)
    plt.show()
    surfpattern(U_record_1D, [x_grid, reduced_t_grid], 'linear',  morphogen=1, rate=0, savefig=False,filename='',logResults=False,normalize=False)
    plt.show()
    # query = simulationOutput_to_sql(simulation_param_dict, model_param_dict,U_final_1D,U_record_1D,ssID=ssID)


#%%
no_growth = False
if no_growth == True:
    growth='nogrowth'
    boundaryCoeff=1
    simulation_param_dict = {'L':L, 'dx':dx, 'J':J, 'T':T, 'dt':dt, 'N':N, 'boundaryCoeff':boundaryCoeff, 'growth':growth, 'growth rate': rate}
    print(simulation_param_dict)

    st = time.time()
    U_final,U_record, U0, x_grid, reduced_t_grid= cn_nogrowth(par_dict,L,J,T,N, circuit_n,boundaryCoeff=boundaryCoeff, tqdm_disable=False)
    elapsed_time = time.time() - st
    print('Execution time:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
    #plot
    plot1D(U_final, savefig=False,filename='')
    plt.show()
    surfpattern(U_record, [x_grid, reduced_t_grid], 'linear',morphogen=0, rate=0, savefig=False,filename='',logResults=False,normalize=False)
    plt.show()
    surfpattern(U_record, [x_grid, reduced_t_grid], 'linear',  morphogen=1, rate=0, savefig=False,filename='',logResults=False,normalize=False)
    plt.show()

    query = simulationOutput_to_sql(simulation_param_dict, model_param_dict,U_final_1D,U_record_1D)



#%%
open_boundary = False
if open_boundary == True:
    growth='nogrowth'
    boundaryCoeff=1
    
    simulation_param_dict = {'L':L, 'dx':dx, 'J':J, 'T':T, 'dt':dt, 'N':N, 'boundaryCoeff':boundaryCoeff, 'growth':growth, 'growth rate': rate}
    
    st = time.time()
    U_final,U_record, U0, x_grid, reduced_t_grid= cn_nogrowth(par_dict,L,J,T,N, circuit_n, boundaryCoeff=2,tqdm_disable=False)
    elapsed_time = time.time() - st
    print('Execution time:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
    #plot
    plot1D(U_final, savefig=False,filename='')
    plt.show()
    surfpattern(U_record, [x_grid, reduced_t_grid], 'linear',morphogen=0, rate=0, savefig=False,filename='',logResults=False,normalize=False)
    plt.show()
    surfpattern(U_record, [x_grid, reduced_t_grid], 'linear',  morphogen=1, rate=0, savefig=False,filename='',logResults=False,normalize=False)
    plt.show()

    query = simulationOutput_to_sql(simulation_param_dict, model_param_dict,U_final_1D,U_record_1D)

# print(np.sum(U_final[0]))


# %%
