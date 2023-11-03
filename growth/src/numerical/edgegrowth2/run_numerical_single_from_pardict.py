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
from numerical.cn_nogrowth import cn_nogrowth


from numerical.cn_plot import plot1D, surfpattern

import pickle
import matplotlib.pyplot as plt
import time
import numpy as np

#system parameters
circuit_n = 'turinghill'
variant=0
n_samples = 1000000

# par_dicxt = {'c1':0.1, 'c2':1,'c3':0.9,'c4':1, 'd_A': 1, 'd_B':10}
# df= pickle.load( open(modellingpath + "/growth/out/analytical/turing/turing_df_%s_variant%r_%rparametersets.pkl"%(circuit_n,variant,n_samples), "rb"))
# df= pickle.load( open(modellingpath + "/growth/out/analytical/turing/turing_df_%s_variant%r_%rparametersets.pkl"%(circuit_n,variant,n_samples), "rb"))

# df = multiple_df.xs(0, level=1)
#solver parameters
L=10; dx =0.1; J = int(L/dx)
T =20; dt = 0.02; N = int(T/dt)


#solver parameters
L=50; dx =0.05; J = int(L/dx)
T =2000; dt = 0.005; N = int(T/dt)
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
parID=544548 ;ssID=0#parameter set to use
# par_dict = df.loc[parID,ssID].to_dict()
# ssID=par_dict['']

# print(par_dict)

# print(par_dict)
# par_dict = df.loc[parID].to_dict()
# print(f'estimated wavelenght = {par_dict["estimated_wvl"]}')


# print(par_dict)
# model_param_dict = {'parID':parID, 'circuit_n':circuit_n,'variant':variant, 'n_samples':n_samples}



#%%
#run
ss_list= np.array([0.02938385, 1.1829765 ])
par_dict= {'ss_list':ss_list,'ba': 0.01, 'bb': 0.01, 'Va': 887.3976810634131, 'Vb': 621.5082654325614, 'kaa': 0.21239643530338062, 'kba': 0.11104906542778901, 'kab': 1.6551470970511883, 'kbb': 0.9597257152142981, 'mua': 5.293891273660762, 'mub': 0.17398259803013874, 'd_B': 0.0014465561532232226, 'd_A': 1.0, 'n': 2.0, 'ss_n': 1.0,  'ss_class': 'unstable point', 'system_class': 'simple unstable', 'maxeig': (3.5570492866984096+0j)}

nogrowth = True
if nogrowth == True:
    growth = 'nogrowth'
    boundaryCoeff=1
    simulation_param_dict = {'L':L, 'dx':dx, 'J':J, 'T':T, 'dt':dt, 'N':N, 'boundaryCoeff':boundaryCoeff, 'growth':growth, 'growth rate': rate}


    U_final_1D,U_record_1D, U0, x_grid, reduced_t_grid= cn_nogrowth(par_dict,L,J,T,N, circuit_n,boundaryCoeff=1, tqdm_disable=False)
    st = time.time()
    plot1D(U_final_1D, savefig=False,filename='')
    plt.show()
    elapsed_time = time.time() - st
    print('Execution time no numba:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
    surfpattern(U_record_1D, L,dx,J,T, 'linear',morphogen=0, rate=0, savefig=False,filename='',logResults=False,normalize=False)
    plt.show()



#%%
open_boundary = False
if open_boundary == True:
    growth='openboundary'
    boundaryCoeff=2
    simulation_param_dict = {'L':L, 'dx':dx, 'J':J, 'T':T, 'dt':dt, 'N':N, 'boundaryCoeff':boundaryCoeff, 'growth':growth, 'growth rate': rate}
    print(simulation_param_dict)


    U_final_1D,U_record_1D, U0, x_grid, reduced_t_grid= cn_nogrowth(par_dict,L,J,T,N, circuit_n,boundaryCoeff=1, tqdm_disable=False)
    st = time.time()
    plot1D(U_final_1D, savefig=False,filename='')
    plt.show()
    elapsed_time = time.time() - st
    print('Execution time no numba:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
    surfpattern(U_record_1D, L,dx,J,T, 'linear',morphogen=0, rate=0, savefig=False,filename='',logResults=False,normalize=False)
    plt.show()


#%%
edgegrowth2 = True
if edgegrowth2 == True:
    growth='edgegrowth2'
    boundaryCoeff=2
    
    simulation_param_dict = {'L':L, 'dx':dx, 'J':J, 'T':T, 'dt':dt, 'N':N, 'boundaryCoeff':boundaryCoeff, 'growth':growth, 'growth rate': rate}

   
    mechanism = 'edgegrowth2'
    boundaryCoeff=2
    st = time.time()
    U_final_1D,U_record_1D, U0, x_grid, reduced_t_grid, cellMatrix= cn_edgegrowth2_numba(par_dict,L,J,T,N, circuit_n, rate=rate, boundaryCoeff=boundaryCoeff, tqdm_disable=False)
    elapsed_time = time.time() - st
    print('Execution time:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
    plot1D(U_final_1D, savefig=False,filename='')
    plt.show()
    elapsed_time = time.time() - st
    print('Execution time no numba:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
    surfpattern(U_record_1D, L,dx,J,T, 'linear',morphogen=0, rate=0, savefig=False,filename='',logResults=False,normalize=False)
    plt.show()

# %%
