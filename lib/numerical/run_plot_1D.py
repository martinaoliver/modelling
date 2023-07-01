
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

def simulate_plot_growth(par_dict,L,J,T,N, circuit_n, rate, tqdm_disable=False):
    st = time.time()
    U_final,U_record, U0, x_grid, reduced_t_grid, cellMatrix= cn_edgegrowth2(par_dict,L,J,T,N, circuit_n, rate=rate, boundaryCoeff=2, tqdm_disable=tqdm_disable)
    elapsed_time = time.time() - st
    print('Execution time no numba:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
    plt.scatter(x_grid,cellMatrix)
    plt.show()
    

    #plot
    plot1D(U_final, savefig=False,filename='')
    plt.show()
    surfpattern(U_record, [x_grid, reduced_t_grid], 'linear',morphogen=0, rate=0, savefig=False,filename='',logResults=False,normalize=False)
    plt.show()
    surfpattern(U_record, [x_grid, reduced_t_grid], 'linear',  morphogen=1, rate=0, savefig=False,filename='',logResults=False,normalize=False)
    plt.show()

    return U_final, U_record


def simulate_plot_nogrowth(par_dict,L,J,T,N, circuit_n, tqdm_disable=False):
    st = time.time()
    U_final,U_record, U0, x_grid, reduced_t_grid= cn_nogrowth(par_dict,L,J,T,N, circuit_n, tqdm_disable=tqdm_disable)
    elapsed_time = time.time() - st
    print('Execution time:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
    #plot
    plot1D(U_final, savefig=False,filename='')
    plt.show()
    surfpattern(U_record, [x_grid, reduced_t_grid], 'linear',morphogen=0, rate=0, savefig=False,filename='',logResults=False,normalize=False)
    plt.show()
    surfpattern(U_record, [x_grid, reduced_t_grid], 'linear',  morphogen=1, rate=0, savefig=False,filename='',logResults=False,normalize=False)
    plt.show()
    return U_final, U_record


def simulate_plot_nogrowth_numba(par_dict,L,J,T,N, circuit_n, tqdm_disable=False):
    st = time.time()
    U_final,U_record, U0, x_grid, reduced_t_grid= cn_nogrowth_numba(par_dict,L,J,T,N, circuit_n, tqdm_disable=tqdm_disable)
    elapsed_time = time.time() - st
    print('Execution time:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
    #plot
    plot1D(U_final, savefig=False,filename='')
    plt.show()
    surfpattern(U_record, [x_grid, reduced_t_grid], 'linear',morphogen=0, rate=0, savefig=False,filename='',logResults=False,normalize=False)
    plt.show()
    surfpattern(U_record, [x_grid, reduced_t_grid], 'linear',  morphogen=1, rate=0, savefig=False,filename='',logResults=False,normalize=False)
    plt.show()
    return U_final, U_record

