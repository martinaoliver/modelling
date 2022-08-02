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

from numerical.cn_edgegrowth import cn_edgegrowth
from numerical.cn_plot import plot1D, surfpattern

#system parameters
circuit_n = 'schnakenberg'
par_dict = {'c1':0.1, 'c2':1,'c3':0.9,'c4':1, 'd_A': 1, 'd_B':10}

#solver parameters
L=100; x_gridpoints=5; J=L*x_gridpoints;I=J 
T=100; t_gridpoints = 20; N=T*t_gridpoints #Number of timepoints

#run
U,U_record, U0, x_grid, reduced_t_grid= cn_edgegrowth(par_dict,L,J,T,N, circuit_n)

#plot
plot1D(U, savefig=False,filename='')
surfpattern(U_record, [x_grid, reduced_t_grid], 'linear',morphogen=1, rate=0, savefig=False,filename='',logResults=False,normalize=False)
surfpattern(U_record, [x_grid, reduced_t_grid], 'linear',  morphogen=0, rate=0, savefig=False,filename='',logResults=False,normalize=False)
