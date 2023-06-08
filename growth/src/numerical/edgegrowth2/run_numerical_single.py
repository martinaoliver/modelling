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
import pickle
import matplotlib.pyplot as plt
import time
import numpy as np

#system parameters
circuit_n = 'turinghill'
variant=9
n_param_sets = 2000000

# par_dicxt = {'c1':0.1, 'c2':1,'c3':0.9,'c4':1, 'd_A': 1, 'd_B':10}
df= pickle.load( open(modellingpath + "/growth/out/analytical/turing/turing_df_%s_variant%r_%rparametersets.pkl"%(circuit_n,variant,n_param_sets), "rb"))

# df = multiple_df.xs(0, level=1)
#solver parameters
L=10; dx =0.1; J = int(L/dx)
T =100; dt = 0.02; N = int(T/dt)

# T =20000; dt = 0.02; N = int(T/dt)
boundaryCoeff=1;rate=L/T

suggesteddt = float(dx*dx*2)

print(f'suggested dt = {suggesteddt}, used dt = {dt}')

suggesteddt = float(dx*dx*2)
print(f'suggested dt = {suggesteddt}, used dt = {dt}')

# parID= (14414,0) #parameter set to use
parID= (3) #parameter set to use
par_dict = df.iloc[parID].to_dict()
print(par_dict)
parameter_to_modify = ['ba', 'bb', 'Va', 'Vb', 'mua', 'mub']
for parameter in parameter_to_modify:
    par_dict[parameter] = par_dict[parameter] 
print(par_dict)
# par_dict = df.loc[parID].to_dict()
print(f'estimated wavelenght = {par_dict["estimated_wvl"]}')


print(par_dict)

#%%
#run

growth = False
if growth == True:
    st = time.time()
    U,U_record, U0, x_grid, reduced_t_grid, cellMatrix= cn_edgegrowth2(par_dict,L,J,T,N, circuit_n, rate=rate, boundaryCoeff=2)
    elapsed_time = time.time() - st
    print('Execution time no numba:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
    plt.scatter(x_grid,cellMatrix)
    plt.show()
    

    #plot
    plot1D(U, savefig=False,filename='')
    plt.show()
    surfpattern(U_record, [x_grid, reduced_t_grid], 'linear',morphogen=0, rate=0, savefig=False,filename='',logResults=False,normalize=False)
    plt.show()
    surfpattern(U_record, [x_grid, reduced_t_grid], 'linear',  morphogen=1, rate=0, savefig=False,filename='',logResults=False,normalize=False)
    plt.show()

#%%
no_growth = True
if no_growth == True:
    st = time.time()
    U_final,U_record, U0, x_grid, reduced_t_grid= cn_nogrowth(par_dict,L,J,T,N, circuit_n, tqdm_disable=False)
    elapsed_time = time.time() - st
    print('Execution time:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
    #plot
    plot1D(U_final, savefig=False,filename='')
    plt.show()
    surfpattern(U_record, [x_grid, reduced_t_grid], 'linear',morphogen=0, rate=0, savefig=False,filename='',logResults=False,normalize=False)
    plt.show()
    surfpattern(U_record, [x_grid, reduced_t_grid], 'linear',  morphogen=1, rate=0, savefig=False,filename='',logResults=False,normalize=False)
    plt.show()


print(np.sum(U_final[0]))


#%%
no_growth_numba = False
if no_growth_numba == True:
    st = time.time()
    U_final,U_record, U0, x_grid, reduced_t_grid= cn_nogrowth_numba(par_dict,L,J,T,N, circuit_n, tqdm_disable=False)
    elapsed_time = time.time() - st
    print('Execution time:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
    #plot
    plot1D(U_final, savefig=False,filename='')
    plt.show()
    surfpattern(U_record, [x_grid, reduced_t_grid], 'linear',morphogen=0, rate=0, savefig=False,filename='',logResults=False,normalize=False)
    plt.show()
    surfpattern(U_record, [x_grid, reduced_t_grid], 'linear',  morphogen=1, rate=0, savefig=False,filename='',logResults=False,normalize=False)
    plt.show()


print(np.sum(U_final[0]))


# %%
