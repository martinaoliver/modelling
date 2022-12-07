#############
###paths#####
#############
import sys
import os

from importlib_metadata import distribution
pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
sys.path.append(modellingpath + '/lib')
#############

from numerical.cn_edgegrowth2_numba import cn_edgegrowth2 as cn_edgegrowth2_numba
from numerical.cn_edgegrowth2 import cn_edgegrowth2
from numerical.cn_nogrowth import cn_nogrowth
from numerical.cn_plot import plot1D, surfpattern
import pickle
import matplotlib.pyplot as plt
import time
import numpy as np
from tqdm import tqdm

#system parameters
circuit_n = 'turinghill'
variant=0 
n_param_sets = 2000000

# par_dict = {'c1':0.1, 'c2':1,'c3':0.9,'c4':1, 'd_A': 1, 'd_B':10}
df= pickle.load( open(modellingpath + "/growth/out/analytical/turing/turing_df_%s_variant%r_%rparametersets.pkl"%(circuit_n,variant,n_param_sets), "rb"))
# df= pickle.load( open(modellingpath + "/growth/out/analytical/lsa_dataframes/lsa_df_%s_variant%r_%rparametersets.pkl"%(circuit_n,variant,n_param_sets), "rb"))
# df= pickle.load( open(modellingpath + "/growth/out/analytical/instability/instability_df_%s_variant%r_%rparametersets.pkl"%(circuit_n,variant,n_param_sets), "rb"))
df.index.names = ['parID','ss']
# df = multiple_df.xs(0, level=1)
#solver parameters

# L=250; x_gridpoints=1; J=L*x_gridpoints;I=J 
# T=15000; t_gridpoints = 5; N=T*t_gridpoints #Number of timepoints <below 3 is bad if x_gridpoints=1
# boundaryCoeff=2;rate=0.01

# L=500; x_gridpoints=2; J=L*x_gridpoints;I=J 
# T=3000; t_gridpoints = 5; N=T*t_gridpoints #Number of timepoints <below 3 is bad if x_gridpoints=1
# boundaryCoeff=2;rate=0.1



# L=500; x_gridpoints=1; J=L*x_gridpoints;I=J 
# T=3000; t_gridpoints = 5; N=T*t_gridpoints #Number of timepoints <below 3 is bad if x_gridpoints=1
# boundaryCoeff=2;rate=0.1


L=500; dx =1; J = int(L/dx)
T =100; dt = 0.05; N = int(T/dt)
boundaryCoeff=2;rate=0.1

L=50; dx =1; J = int(L/dx)
T =500; dt = 0.005; N = int(T/dt)
boundaryCoeff=2;rate=0.1


print(len(df))

filename= lambda mechanism, parID: 'circuit%s_variant%s_bc%s_%s_rate%s_ID%s_L%r_J%r_T%r_N%r'%(circuit_n,variant,boundaryCoeff, mechanism,rate,parID,L,J,T,N)


for parID,ss in tqdm(df.index[:1], disable=False):
# for parID,ss in tqdm([(1271426,4)], disable=True):
    parIDss = f'{parID}.{ss}'
    # print(parIDss)
    mechanism = 'edgegrowth2'
    par_dict = df.loc[(parID,ss)].to_dict()
    # print(par_dict)
    # print(mechanism)
    U_final,U_record, U0, x_grid, reduced_t_grid, cellMatrix= cn_edgegrowth2_numba(par_dict,L,J,T,N, circuit_n, rate=rate, boundaryCoeff=boundaryCoeff, tqdm_disable=False)
    # plt.scatter(x_grid,cellMatrix)
    # plt.show()
    surfpattern(U_record, [x_grid, reduced_t_grid], 'linear',morphogen=0, rate=0, savefig=False,filename='',logResults=False,normalize=False)
    plt.show()
    # pickle.dump(U_final, open(modellingpath + '/growth/out/numerical/%s/%s/simulation/2Dfinal_%s.pkl'%(circuit_n,mechanism,filename(mechanism,parIDss)), 'wb'))
    # pickle.dump(U_record, open(modellingpath + '/growth/out/numerical/%s/%s/simulation/2Drecord_%s.pkl'%(circuit_n,mechanism,filename(mechanism,parIDss)), 'wb'))
    plot1D(U_final)
    plt.show()
    mechanism = 'nogrowth'
    # print(mechanism)

    U_final,U_record, U0, x_grid, reduced_t_grid= cn_nogrowth(par_dict,L,J,T,N, circuit_n, tqdm_disable=False)
    # pickle.dump(U_final, open(modellingpath + '/growth/out/numerical/%s/%s/simulation/2Dfinal_%s.pkl'%(circuit_n,mechanism,filename(mechanism,parIDss)), 'wb'))
    # pickle.dump(U_record, open(modellingpath + '/growth/out/numerical/%s/%s/simulation/2Drecord_%s.pkl'%(circuit_n,mechanism,filename(mechanism,parIDss)), 'wb'))
    surfpattern(U_record, [x_grid, reduced_t_grid], 'linear',morphogen=0, rate=0, savefig=False,filename='',logResults=False,normalize=False)
    plt.show()
    plot1D(U_final)
    plt.show()

# #run

# no_numba = False
# if no_numba == True:
#     st = time.time()
#     U,U_record, U0, x_grid, reduced_t_grid, cellMatrix= cn_edgegrowth2(par_dict,L,J,T,N, circuit_n, rate=0.1, boundaryCoeff=2)
#     elapsed_time = time.time() - st
#     print('Execution time numba:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
#     plt.scatter(x_grid,cellMatrix)
#     plt.show()

#     #plot
#     plot1D(U, savefig=False,filename='')
#     plt.show()
#     surfpattern(U_record, [x_grid, reduced_t_grid], 'linear',morphogen=0, rate=0, savefig=False,filename='',logResults=False,normalize=False)
#     plt.show()
#     surfpattern(U_record, [x_grid, reduced_t_grid], 'linear',  morphogen=1, rate=0, savefig=False,filename='',logResults=False,normalize=False)
#     plt.show()


# numba=False
# if numba == True: 
#     st = time.time()
#     U,U_record, U0, x_grid, reduced_t_grid, cellMatrix= cn_edgegrowth2_numba(par_dict,L,J,T,N, circuit_n, rate=0.1, boundaryCoeff=2)
#     elapsed_time = time.time() - st
#     print('Execution time:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))


# no_growth = True
# if no_growth == True:
#     st = time.time()
#     U,U_record, U0, x_grid, reduced_t_grid= cn_nogrowth(par_dict,L,J,T,N, circuit_n, tqdm_disable=False)
#     elapsed_time = time.time() - st
#     print('Execution time:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
#     #plot
#     plot1D(U, savefig=False,filename='')
#     plt.show()
#     surfpattern(U_record, [x_grid, reduced_t_grid], 'linear',morphogen=0, rate=0, savefig=False,filename='',logResults=False,normalize=False)
#     plt.show()
#     surfpattern(U_record, [x_grid, reduced_t_grid], 'linear',  morphogen=1, rate=0, savefig=False,filename='',logResults=False,normalize=False)
#     plt.show()


# print(np.sum(U_final[0]))