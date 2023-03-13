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
# import matplotlib.pyplot as plt
# from matplotlib import colors
# from mpl_toolkits.axes_grid1 import make_axes_locatable


#############################
# Specify name of circuit and variant investigated
# circuit_n=14;variant='fitted1';n_species=6
circuit_n=14;variant='2nd';n_species=6
# Specifiy number of parameter sets in parameterset file to be loaded
# balance = 'balanced'
# folder = 'circuit14variantfitted1'
folder ='circuit14variant2ndBalancedTuring'
modelArgs = [circuit_n,variant,n_species,folder]

# Specifiy number of parameter sets in parameterset file to be loaded

# specify dimensions of system
# specify dimensions of system
L=9; dx =0.05; J = int(L/dx)
T =50; dt = 0.05; N = int(T/dt)
boundarycoeff = 1
divisionTimeHours=0.5
p_division=1;seed=1

shape = 'ca'

x_gridpoints=int(1/dx)

# filename= lambda parID: 'circuit%r_variant%s_bc%s_%s_ID%s_L%r_J%r_T%r_N%r'%(circuit_n,variant,boundarycoeff, shape,parID,L,J,T,N)
degDiv=4;filename= lambda parID: 'circuit%r_variant%s_bc%s_%s_ID%s_L%r_J%r_T%r_N%r_degDiv%s'%(circuit_n,variant,boundarycoeff, shape,parID,L,J,T,N,degDiv)
data_path = modellingephemeral + '/3954/paper/out/numerical/colonies/simulation/%s'%(folder)
parID_list = pickle.load( open(data_path + '/parID_list_%s.pkl'%(filename('x')), "rb" ) )

start = 0
stop = len(parID_list) 
# stop=10
plotAllFunctionColonies(parID_list, circuit_n, shape, filename, L,x_gridpoints,folder=folder,start=start, stop=stop, tqdm_disable=False, saveFig=True)

