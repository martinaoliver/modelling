print('hehe')
#############
###paths#####
#############
import sys
import os

pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
modellingephemeral = '/rds/general/ephemeral/user/mo2016/ephemeral/Documents/modelling'
sys.path.append(modellingpath + '/lib')
#############

from numerical.plotAllFunctions import plotAllFunctionColonies

import pickle
import numpy as np
print('hehe2')
# import matplotlib.pyplot as plt
# from matplotlib import colors
# from mpl_toolkits.axes_grid1 import make_axes_locatable


#############################
# Specify name of circuit and variant investigated
# circuit_n=14;variant='fitted1';n_species=6
# circuit_n=14;variant='2nd';n_species=6
# circuit_n=14;variant=int(sys.argv[1]);n_species=6;nsr=0.05
# circuit_n=14;variant=195238;n_species=6;nsr=0.05
circuit_n=14;variant='2nd';n_species=6; Kce=100
# circuit_n=14;variant='fitted7_gaussian4187715_nsr0.01';n_species=6
n_samples=1000000
# n_samples=2000
circuit_n=14;variant='2nd';n_species=6; Kce=100
# circuit_n=14;variant='fitted7_gaussian4187715_nsr0.01';n_species=6
n_samples=1000000
# n_samples=2000
# Specifiy number of parameter sets in parameterset file to be loaded
# balance = 'balanced'
# folder = 'circuit14variantfitted1'
# folder ='circuit14variant2ndBalancedTuring'
# folder = f'circuit14variant{variant}'
# folder = 'circuit14variant2ndBalancedKce100'
folder = 'circuit14variant2ndBalancedKce100'
# folder = 'circuit14variantfitted7_gaussian4187715'
ssID=0
model_param_dict = {'circuit_n':circuit_n,'variant':variant, 'n_samples':n_samples}

modelArgs = [circuit_n,variant,n_species,folder]

# Specifiy number of parameter sets in parameterset file to be loaded


# specify dimensions of system

# # slow
# L=20; dx =0.1; J = int(L/dx)
# T =100; dt = 0.02; N = int(T/dt)
# boundaryCoeff = 2
# division_time_hours=0.5
# p_division=0.38;seed=1
# # slow
# L=20; dx =0.1; J = int(L/dx)
# T =100; dt = 0.02; N = int(T/dt)
# boundaryCoeff = 2
# division_time_hours=0.5
# p_division=0.38;seed=1


# # medium
# L=20; dx =0.1; J = int(L/dx)
# T =50; dt = 0.02; N = int(T/dt)
# boundaryCoeff = 1
# division_time_hours=0.5
# p_division=1;seed=1



# fast
L=20; dx =0.1; J = int(L/dx)
T =25; dt = 0.02; N = int(T/dt)
boundaryCoeff = 1
division_time_hours=0.2
p_division=0.7;seed=1
# fast
L=20; dx =0.1; J = int(L/dx)
T =25; dt = 0.02; N = int(T/dt)
boundaryCoeff = 1
division_time_hours=0.2
p_division=0.7;seed=1

shape = 'ca'

x_gridpoints=int(1/dx)


# filename= lambda parID: 'circuit%r_variant%s_bc%s_%s_ID%s_L%r_J%r_T%r_N%r'%(circuit_n,variant,boundarycoeff, shape,parID,L,J,T,N)
# filename= lambda parID: 'circuit%r_variant%snsr%s_bc%s_%s_ID%s_L%r_J%r_T%r_N%r'%(circuit_n,variant,nsr,boundarycoeff, shape,parID,L,J,T,N)
# filename= lambda parID: 'circuit%r_variant%s_%sparametersets_balanced_Kce%s_bc%s_%s_ID%s_L%r_J%r_T%r_N%r'%(circuit_n,variant,n_samples,Kce,boundarycoeff, shape,parID,L,J,T,N)
data_path = modellingpath + '/3954/paper/out/numerical/colonies/simulation/%s'%(folder)
parID_list = pickle.load( open(data_path + '/parID_list_%s.pkl'%(filename('x')), "rb" ) )

start=0
# start = 825

stop = len(parID_list) 
# stop = 10

plotAllFunctionColonies(parID_list, circuit_n, shape, filename, L,x_gridpoints,folder=folder,start=start, stop=stop, tqdm_disable=False, saveFig=True)

plotAllFunctionColonies_differentSnapshot(parID_list, circuit_n, shape,snapshot, filename, L,x_gridpoints,start=0, stop=10,folder=None, modellingpath=modellingpath, saveFig=True,dpi=2000, tqdm_disable=True, print_parID=False):
